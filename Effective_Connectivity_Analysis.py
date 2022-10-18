import networkx
import pandas as pd
import numpy as np
from scipy.stats import sem
import seaborn as sns
import matplotlib.pyplot as plt
import read_connectivity
import networkx as Nx
from bokeh.plotting import from_networkx
from bokeh.models import Arrow, VeeHead
from bokeh.models import Range1d, Circle, ColumnDataSource, MultiLine
from bokeh.plotting import figure
from networkx.algorithms import community

################# Constants #################
node_regions = [8, 5, 8, 8, 7, 0, 3, 3, 9, 9, 10, 10, 10, 7, 3, 3, 8, 8, 0, 8, 3, 7, 5, 3, 9, 9, 10, 10, 6, 7, 3, 3, 1,
                2, 1, 1, 1, 1, 7, 1, 6, 0, 0, 3, 3, 1, 4, 1, 1, 2]
region_names = np.array(['PMd (6DR)', 'PMd', 'SMA', 'preSMA', 'dlPFC (46D)', 'FEF (8AD)', 'WM', 'ACC', 'Cd',
                         'OFC', 'vmPFC'])
spectral11 = np.array(['#9e0142', '#d53e4f', '#f46d43', '#fdae61', '#fee08b', '#ffffbf', '#e6f598', '#abdda4',
                       '#66c2a5', '#3288bd', '#5e4fa2'])
################# Functions #################


def weights_to_edge_df(A):
    """Transforms a list of nodes and array of parameters into an edge dataframe.
        Inputs:
          A (2D array) - connectivity weights from DCM
        Outputs:
          df (data frame) - edges defined by {Source, Target, Weight}
        """
    nodes = list(range(1, A.shape[0]+1))
    start_nodes = [n for n in nodes for _ in range(len(nodes))]
    end_nodes = nodes * len(nodes)
    A = np.ndarray.flatten(A)
    dict = {'Source': start_nodes, 'Target': end_nodes, 'Weight': A}
    df = pd.DataFrame(data=dict)
    return df


def clean_df(df, threshold=None):
    """Removes edges with a weight lower than a specified threshold.
        Inputs:
          df (data frame) - edge defined by {Source, Target, Weight}
          threshold (number) - minimum edge weight
        Outputs:
          df (data frame) - edges where Weight > threshold
        """
    return df.loc[abs(df['Weight']) > threshold]


def arrow_renderer(df, N, G, plot, edge_sign=True, edge_magnitude=True, self_arrows=True, node_rad=0.25):
    """Draws directed arrows as edges on a plot. Colors them by edge weight sign and sets line width by edge weight
    magnitude if desired.
        Inputs:
          df (data frame) - edges defined by {Source, Target, Weight}
          N (network) - networkx network used to retrieve nodes
          G (graph) - bokeh graph object used to retrieve node coordinates
          plot (figure) - figure to plot the arrows on
          edge_sign (bool) - True if coloring by edge weight sign is desired
          edge_magnitude (bool) - True if line width by edge weight magnitude is desired
          self_arrows (bool) - True if circular arrows for self-affecting nodes are desired
          node_rad (num) - radius of node glyphs in units of axes
        Outputs:
          plot (figure) - figure with arrows added
        """
    coords = G.layout_provider.graph_layout         # key - nodes: value - [x,y]
    # Map node (x,y) to their occurrences as source and target
    df['Source_x'] = [coords[k][0] for k in df['Source']]
    df['Source_y'] = [coords[k][1] for k in df['Source']]
    df['Target_x'] = [coords[k][0] for k in df['Target']]
    df['Target_y'] = [coords[k][1] for k in df['Target']]
    dfself = df[(df['Source']-df['Target'] == 0)]                   # subset of edges that are self-affecting
    df = df[(df['Source']-df['Target'] != 0)]                       # subset of edges that are not self-affecting
    df.reset_index(drop=True, inplace=True)
    max_width = 2
    min_width = 0.5
    # Make arrow head of straight edges terminate on node edge
    n_shift = node_rad / 2 * (2 ** 0.5)
    for i in range(len(df)):                                        # for every edge
        if df.loc[i, 'Target_x'] > df.loc[i, 'Source_x']:           # if arrow points right
            df.loc[i, 'Target_x'] -= n_shift                        # move arrow head left by node_rad
        else:                                                       # if arrow points left
            df.loc[i, 'Target_x'] += n_shift                        # move arrow head right by node_rad
        if df.loc[i, 'Target_y'] > df.loc[i, 'Source_y']:           # if arrow points up
            df.loc[i, 'Target_y'] -= n_shift                        # move arrow head down by node_rad
        else:                                                       # if arrow points down
            df.loc[i, 'Target_y'] += n_shift                        # move arrow head up by node_rad
    # Data source for straight arrows
    source = ColumnDataSource(
        {'Source_x': df['Source_x'], 'Source_y': df['Source_y'], 'Target_x': df['Target_x'], 'Target_y': df['Target_y']})
    arrows = Arrow(end=VeeHead(size=8), source=source, x_start='Source_x', y_start='Source_y', x_end='Target_x',
                   y_end='Target_y', level='underlay')
    if edge_sign:
        df['Color'] = df['Weight'].copy()
        df['Color'].mask(df['Color'] < 0, 'Red', inplace=True)       # red for negative weights
        df['Color'].mask(df['Color'] != 'Red', 'Blue', inplace=True) # blue for positive weights
        source.add(data=df['Color'], name='Color')                   # add color attribute to straight arrow data source
        arrows.line_color = 'Color'                                  # color arrow heads and lines accordingly
        arrows.end.line_color = 'Color'
        arrows.end.fill_color = 'Color'
    if edge_magnitude:
        # Map weights onto new range specified by min and max width
        df['Width'] = (abs(df['Weight']) - min(abs(df['Weight'])))/(max(df['Weight']) - min(abs(df['Weight']))) * \
                      (max_width-min_width) + min_width
        source.add(data=df['Width'], name='Width')                   # add width attribute to straight arrow data source
        arrows.line_width = 'Width'                                  # set arrow line width accordingly
    if self_arrows and dfself.size != 0:
        dfself.reset_index(drop=True, inplace=True)
        circle_rad = 0.9                                         # radius of circular arrow in axes units
        c_shift = circle_rad/2*(2**0.5)                             # circle edge goes through center of node
        ah_offset = node_rad/3                                      # arrow head shift to go through circle
        ah_length = 0.1                                             # set arrow direction
        shift = np.array([-c_shift, c_shift])                       # circle center shift options
        # [[[right to node edge, right to node edge], [left to node edge, left to node edge]],
        # [[right to node center + node_rad/2, right to node center + node_rad/2 ],
        # [left to node center + node_rad/2, left to node center + node_rad/2]]]
        end_x = np.array([[[c_shift-node_rad, c_shift-node_rad], [-(c_shift-node_rad), -(c_shift-node_rad)]],
                          [[c_shift+ah_offset, c_shift+ah_offset], [-c_shift-ah_offset, -c_shift-ah_offset]]])
        # [[[up to node center + node_rad/2, down to node center + node_rad/2],
        # [up to node center + node_rad/2, down to node center + node_rad/2]],
        # [[up to node edge, down to node edge], [up to node edge, down to node edge]]]
        end_y = np.array([[[c_shift+ah_offset, -(c_shift+ah_offset)], [c_shift+ah_offset, -(c_shift+ah_offset)]],
                          [[(c_shift-node_rad), -(c_shift-node_rad)], [(c_shift-node_rad), -(c_shift-node_rad)]]])
        # [[[0.1 left from node edge, 0.1 left from node edge], [0.1 right from node edge, 0.1 right from node edge]],
        # [[vertical, vertical], [vertical, vertical]]]
        start_x = np.array([[[c_shift-node_rad-ah_length, c_shift-node_rad-ah_length],
                             [-(c_shift-node_rad)+ah_length, -(c_shift-node_rad)+ah_length]],
                            [[c_shift+ah_offset, c_shift+ah_offset], [-c_shift-ah_offset, -c_shift-ah_offset]]])
        # [[[horizontal, horizontal], [horizontal, horizontal]],
        # [[0.1 down from node edge, 0.1 up from node edge], [0.1 down from node edge, 0.1 up from node edge]]]
        start_y = np.array([[[c_shift+ah_offset, -(c_shift+ah_offset)], [c_shift+ah_offset, -(c_shift+ah_offset)]],
                            [[(c_shift-node_rad-ah_offset), -(c_shift-node_rad)+ah_offset],
                             [(c_shift-node_rad-ah_offset), -(c_shift-node_rad)+ah_offset]]])
        # Set circle and arrow head locations based on node quadrant
        dfself['x_half'] = (dfself['Source_x'] > 0).astype('uint8')                 # 0 for x<=0 and 1 for x>0
        dfself['y_half'] = (dfself['Source_y'] > 0).astype('uint8')                 # 0 for y<=0 and 1 for y>0
        dfself['angle'] = (abs(dfself['Source_y']) > abs(dfself['Source_x'])).astype('uint8')
        dfself['Source_x'] += shift[dfself['x_half']]                               # shift right for pos x, left for neg
        dfself['Source_y'] += shift[dfself['y_half']]                               # shift up for pos y, down for neg
        # Determine arrow head locations as shifts from circle center using arrays and node quadrant/angle
        dfself['ar_end_x'] = dfself['Source_x'] + end_x[dfself['angle'], dfself['x_half'], dfself['y_half']]
        dfself['ar_end_y'] = dfself['Source_y'] + end_y[dfself['angle'], dfself['x_half'], dfself['y_half']]
        dfself['ar_start_x'] = dfself['Source_x'] + start_x[dfself['angle'], dfself['x_half'], dfself['y_half']]
        dfself['ar_start_y'] = dfself['Source_y'] + start_y[dfself['angle'], dfself['x_half'], dfself['y_half']]
        # Data source for circular arrows
        circle_source = ColumnDataSource({'Source_x': dfself['Source_x'], 'Source_y': dfself['Source_y'],
                                          'ar_start_x': dfself['ar_start_x'], 'ar_start_y': dfself['ar_start_y'],
                                          'ar_end_x': dfself['ar_end_x'], 'ar_end_y': dfself['ar_end_y']})
        # initialize circles and arrow heads objects
        circles = plot.circle(x='Source_x', y='Source_y', source=circle_source, fill_color=None, radius=circle_rad, level='underlay')
        arrow_heads = Arrow(end=VeeHead(size=8), source=circle_source,
                            x_start='ar_start_x', y_start='ar_start_y',x_end='ar_end_x', y_end='ar_end_y',
                            level='underlay')
        if edge_sign:                                               # same as straight arrows
            dfself['Color'] = dfself['Weight']
            dfself['Color'].mask(dfself['Color'] < 0, 'Red', inplace=True)
            dfself['Color'].mask(dfself['Color'] != 'Red', 'Blue', inplace=True)
            circle_source.add(data=dfself['Color'], name='Color')   # add color to circular arrow data source
            arrow_heads.end.line_color = 'Color'                    # color arrow heads accordingly
            arrow_heads.end.fill_color = 'Color'
            # can't update attributes of circles like you can for arrow heads so you have to redefine with line color
            circles = plot.circle(x='Source_x', y='Source_y', source=circle_source, fill_color=None, radius=circle_rad,
                                  level='underlay', line_color='Color')
        if edge_magnitude:                                          # same as straight arrows
            dfself['Width'] = (abs(dfself['Weight']) - min(abs(dfself['Weight']))) / (
                    max(dfself['Weight']) - min(abs(dfself['Weight']))) * (max_width - min_width) + min_width
            circle_source.add(data=dfself['Width'], name='Width')   # add line width to circular arrow data source
            if edge_sign:                                           # if edge color and weight are desired
                # redefine circles with color and width attributes
                circles = plot.circle(x='Source_x', y='Source_y', source=circle_source, fill_color=None, radius=circle_rad,
                                      level='underlay', line_color='Color', line_width='Width')
            else:                                                   # if only edge weight is desired
                # redefine circles with only line width attributes
                circles = plot.circle(x='Source_x', y='Source_y', source=circle_source, fill_color=None, radius=circle_rad,
                                      level='underlay', line_width='Width')
        plot.add_layout(circles)
        plot.add_layout(arrow_heads)
    plot.add_layout(arrows)
    return plot


def signed_in_degree(df):
    """ Computes the excitatory (+) and inhibitory (-) weighted in degrees for each node.
        Inputs:
          df (data frame) - edges defined by {Source, Target, Weight}
        Outputs:
          exc_in_degree (dict) - node: weighted in degree of positive edges
          inh_in_degree (dict) - node: weighted in degree of negative edges
        """
    excitatory_df = df.copy()                                                     # copy so df doesn't change
    excitatory_df['Weight'].mask(excitatory_df['Weight'] < 0, 0, inplace=True)    # set negative weights to 0
    excitatory_N = Nx.from_pandas_edgelist(excitatory_df, 'Source', 'Target', 'Weight', create_using=Nx.DiGraph())
    exc_in_degree = dict(excitatory_N.in_degree(weight='Weight'))
    inhibitory_df = df.copy()
    inhibitory_df['Weight'].mask(inhibitory_df['Weight'] > 0, 0, inplace=True)    # set positive weights to 0
    inhibitory_N = Nx.from_pandas_edgelist(inhibitory_df, 'Source', 'Target', 'Weight', create_using=Nx.DiGraph())
    inh_in_degree = dict(inhibitory_N.in_degree(weight='Weight'))
    return exc_in_degree, inh_in_degree


def signed_out_degree(df):
    """ Computes the excitatory (+) and inhibitory (-) weighted out degrees for each node.
        Inputs:
          df (data frame) - edges defined by {Source, Target, Weight}
        Outputs:
          exc_out_degree (dict) - node: weighted out degree of positive edges
          inh_out_degree (dict) - node: weighted out degree of negative edges
        """
    excitatory_df = df.copy()                                                       # copy so df doesn't change
    excitatory_df['Weight'].mask(excitatory_df['Weight'] < 0, 0, inplace=True)      # set negative weights to 0
    excitatory_N = Nx.from_pandas_edgelist(excitatory_df, 'Source', 'Target', 'Weight', create_using=Nx.DiGraph())
    exc_out_degree = dict(excitatory_N.out_degree(weight='Weight'))
    inhibitory_df = df.copy()
    inhibitory_df['Weight'].mask(inhibitory_df['Weight'] > 0, 0, inplace=True)      # set positive weights to 0
    inhibitory_N = Nx.from_pandas_edgelist(inhibitory_df, 'Source', 'Target', 'Weight', create_using=Nx.DiGraph())
    inh_out_degree = dict(inhibitory_N.out_degree(weight='Weight'))
    return exc_out_degree, inh_out_degree


def connectivity_matrix(A, trial, cond, node_regions):
    """ Generates a heatmap connectivity matrix for a parameter matrix A
        Inputs:
          A (2D matrix) - connectivity weights between source and target nodes from DCM
          trial (num) - trial number for plot title
          cond (string) - condition type for BMI task
        """
    nodes = list(range(1, 51))
    sq = np.column_stack((node_regions, nodes, A))                              # stack node regions and labels before A
    sq = sq[sq[:, 0].argsort()]                                                 # sort rows by region
    node_labels = [int(x) for x in sq[:, 1]]                                    # new node label order
    sq = sq[:, 2:]                                                              # retrieve half sorted A
    sq = np.row_stack((sq, nodes, node_regions))                                # stack node labels and regions under A
    sq = sq[sq[-1, :].argsort()]                                                # sort columns by region
    grid_kws = {"height_ratios": (.9, .05), "hspace": .2}                       # format heatmap
    f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws, figsize=[8.4, 9.4])
    ax = sns.heatmap(sq, linewidths=0.5, ax=ax, cmap='RdBu',cbar_ax=cbar_ax, cbar_kws={"orientation": "horizontal"},
                     xticklabels=node_labels, yticklabels=node_labels)
    ax.set_xlabel("Source Node")
    ax.set_ylabel("Target Node")
    ax.set_title(f'Connectivity Weights - Trial {trial} - Condition {cond}')
    plt.show()
    return


def get_circle_layout(N, scale=20, center=(0, 0)):
    """ Generates a dictionary of nodes and coordinates with the nodes segregated by region
        Inputs:
          N (network) - networkx network
          scale (num) - radius of the circle
          center (tuple) - center of the circle
        Outputs:
          circle_layout (dictionary) - key - sorted nodes: value - [x,y]
        """
    dummy_graph = from_networkx(N, Nx.circular_layout, scale=scale, center=center)
    coords = dummy_graph.layout_provider.graph_layout                               # key - nodes: value - [x,y]
    nodes = list(range(1, 51))
    sorted_nodes = [x for y, x in sorted(zip(node_regions, nodes))]                 # nodes sorted by region
    circle_layout = dict(zip(sorted_nodes, coords.values()))                        # key - sorted nodes: value - [x,y]
    return circle_layout


def degree_comparison(trials, condition_split=False, region_group=True):
    CM = read_connectivity.get_CM()
    f = 1                                                           # figure iterator
    cond_map = {'L': 1, 'M': 2, 'H': 3}                             # subplot selector
    if condition_split:
        # key - condition: value - array for degree of each node (50) for any amount of trials
        IE = dict(zip(['L', 'M', 'H'], np.zeros((3, 1, 50))))
        OE = dict(zip(['L', 'M', 'H'], np.zeros((3, 1, 50))))
        II = dict(zip(['L', 'M', 'H'], np.zeros((3, 1, 50))))
        OI = dict(zip(['L', 'M', 'H'], np.zeros((3, 1, 50))))
        trs =dict(zip(['L', 'M', 'H'], [[0], [0], [0]]))           # key - condition: value - list of trials for each
        for t in trials:
            cond = CM['Target'][t]
            df = weights_to_edge_df(CM['A'][t])
            df = clean_df(df, 0.25)
            df.reset_index(drop=True, inplace=True)
            N = Nx.from_pandas_edgelist(df, 'Source', 'Target', 'Weight',
                                        create_using=Nx.DiGraph())          # network - use to get number of nodes
            in_exc, in_inh = signed_in_degree(df)
            out_exc, out_inh = signed_out_degree(df)
            if len(N.nodes) == 50:                                        # if all nodes are present given the threshold
                ie = dict(sorted(in_exc.items()))                         # sort each degree dictionary by node
                ii = dict(sorted(in_inh.items()))
                oe = dict(sorted(out_exc.items()))
                oi = dict(sorted(out_inh.items()))
                # stack the node degrees for trial t in the key=condition value of the degree dictionaries
                IE[cond] = np.row_stack((IE[cond], list(ie.values())))
                II[cond] = np.row_stack((II[cond], list(ii.values())))
                OE[cond] = np.row_stack((OE[cond], list(oe.values())))
                OI[cond] = np.row_stack((OI[cond], list(oi.values())))
                trs[cond].append(trs[cond][-1]+1)                     # running count of the number of trials for cond
        trs = {k: trs[k][1:] for k in trs}                            # pop the zero off the front of each trial list
        matrices = [IE, II, OE, OI]
        measures = ['In Excitatory', 'In Inhibitory', 'Out Excitatory', 'Out Inhibitory']
        for mat, meas in zip(matrices, measures):                           # for each of the four degree measurements
            fig = plt.figure(figsize=[11, 7])                               # make a figure
            plt.xlabel('Trial', labelpad=20)                                # common xlabel
            plt.ylabel(f'{meas} Degree', labelpad=35)                       # common ylabel
            plt.title(f'Shift in Weighted {meas} Degree', pad=25)           # common title
            plt.xticks([])                                                  # suppress common ticks
            plt.yticks([])
            for condition, cond in cond_map.items():                        # for each condition
                ax = fig.add_subplot(3, 1, cond)                            # add a subplot
                ax.title.set_text(f'Condition {condition}')                 # subplot title
                ax.set_xticks(trs[condition])                               # subplot ticks are trials for the condition
                if region_group:                                            # if you want to group nodes by region
                    # key - region code: value - array for degree of nodes in the region for any amount of trials
                    grouped_nodes = dict(zip(list(range(11)), np.zeros((11, len(trs[condition]), 1))))
                    for c in range(50):                                     # for each node
                        # add it to the value of the correct region based on the list of node_regions
                        grouped_nodes[node_regions[c]] = np.column_stack((grouped_nodes[node_regions[c]], mat[condition][1:, c]))
                    grouped_nodes = {r: grouped_nodes[r][:, 1:] for r in grouped_nodes}                    # remove the first column of zeros
                    grouped_nodes_mean = {r: np.mean(grouped_nodes[r], axis=1) for r in grouped_nodes}     # average degree over trials for each region
                    grouped_nodes_sd = {r: sem(grouped_nodes[r], axis=1) for r in grouped_nodes}           # standard error of the mean
                    for r in range(11):                                                                    # for each region
                        plt.plot(trs[condition], grouped_nodes_mean[r], color=spectral11[r])               # plot the degree for applicable number of trials
                        plt.fill_between(trs[condition], grouped_nodes_mean[r]-grouped_nodes_sd[r],
                                         grouped_nodes_mean[r]+grouped_nodes_sd[r], alpha=0.3,
                                         facecolor=spectral11[r])  # shade around line with SEM; coloring is same as network
                else:                                                           # if you don't want to group by region
                    for c in range(50):                                         # for each node
                        plt.plot(trs[condition], mat[condition][1:, c])         # plot degree for applicable number of trials
            axes = fig.get_axes()                                               # force y axes to have same limits
            axes[0].get_shared_y_axes().join(*axes)
            plt.tight_layout()
            plt.subplots_adjust(hspace=0.5, top=0.85)                       # space for titles
            plt.show()                                                      # show plot for one degree measure
    else:                                                                   # w/o considering condition
        IE = np.zeros((1, 50))                                              # frame for degrees of 50 nodes
        OE = np.zeros((1, 50))
        II = np.zeros((1, 50))
        OI = np.zeros((1, 50))
        trs = list()                                                        # empty list for trials with all nodes
        for t in trials:                                                    # for each trial
            df = weights_to_edge_df(CM['A'][t])
            df = clean_df(df, 0.25)
            df.reset_index(drop=True, inplace=True)
            N = Nx.from_pandas_edgelist(df, 'Source', 'Target', 'Weight',
                                        create_using=Nx.DiGraph())
            in_exc, in_inh = signed_in_degree(df)
            out_exc, out_inh = signed_out_degree(df)
            if len(N.nodes) == 50:                                        # if all nodes are present given the threshold
                ie = dict(sorted(in_exc.items()))                         # sort each degree dictionary by node
                ii = dict(sorted(in_inh.items()))
                oe = dict(sorted(out_exc.items()))
                oi = dict(sorted(out_inh.items()))
                # stack the node degrees for trial t in the degree dictionaries
                IE = np.row_stack((IE, list(ie.values())))
                II = np.row_stack((II, list(ii.values())))
                OE = np.row_stack((OE, list(oe.values())))
                OI = np.row_stack((OI, list(oi.values())))
                trs.append(t)                                               # add the trial to the list
        IE = IE[1:, :]                                                      # remove the first row of zeros
        II = II[1:, :]
        OE = OE[1:, :]
        OI = OI[1:, :]
        matrices = [IE, II, OE, OI]
        measures = ['In Excitatory', 'In Inhibitory', 'Out Excitatory', 'Out Inhibitory']
        cond_color = {'L': 'red', 'M': 'yellow', 'H': 'green'}
        for mat, meas in zip(matrices, measures):                           # for each of the four degree measurements
            plt.figure(f)                                                   # make a new figure
            f += 1
            if region_group:  # if you want to group nodes by region
                # key - region code: value - array for degree of nodes in the region for any amount of trials
                grouped_nodes = dict(zip(list(range(11)), np.zeros((11, len(trs), 1))))
                for c in range(50):  # for each node
                    # add it to the value of the correct region based on the list of node_regions
                    grouped_nodes[node_regions[c]] = np.column_stack((grouped_nodes[node_regions[c]], mat[:, c]))
                grouped_nodes = {r: grouped_nodes[r][:, 1:] for r in grouped_nodes}  # remove the first column of zeros
                grouped_nodes_mean = {r: np.mean(grouped_nodes[r], axis=1) for r in
                                      grouped_nodes}  # average degree over trials for each region
                grouped_nodes_sd = {r: sem(grouped_nodes[r], axis=1) for r in
                                    grouped_nodes}                               # standard error of the mean
                for r in range(11):                                              # for each region
                    plt.plot(trs, grouped_nodes_mean[r],
                             color=spectral11[r])  # plot the degree for applicable number of trials
                    plt.fill_between(trs, grouped_nodes_mean[r] - grouped_nodes_sd[r],
                                     grouped_nodes_mean[r] + grouped_nodes_sd[r], alpha=0.3,
                                     facecolor=spectral11[r])  # shade around line with SEM; coloring is same as network
            else:
                for c in range(50):                                              # for each node
                    plt.plot(trs, mat[:, c], color=spectral11[node_regions[c]])  # plot the degree over all trials
            for t in trs:                                                        # for each trial
                plt.axvspan(t-0.5, t+0.5, alpha=0.15, color=cond_color[CM['Target'][t]],
                            linewidth=0)   # vertical shading by condition
            plt.title(f'Shift in Weighted {meas} Degree Over {len(trs)} Trials')
            plt.xlabel('Trial')
            plt.ylabel(f'{meas} Degree')
            plt.show()


def efficiency_comparison(trials, condition_split=False):
    """ Plots the local and global efficiency over trials with w/ or w/o splitting by condition
            Inputs:
              trials (list) - the trials that should be considered
              condition_split (bool) - split the plots by condition or not
            """
    CM = read_connectivity.get_CM()
    cond_map = {'L': 1, 'M': 2, 'H': 3}                                    # subplot selector
    cond_color = {'L': 'red', 'M': 'yellow', 'H': 'green'}                 # condition color map
    if condition_split:
        # key - condition: value - efficiency for any amount of trials
        # GE = Global Efficiency, LE = Local Efficiency
        GE = dict(zip(['L', 'M', 'H'], [[0], [0], [0]]))
        LE = dict(zip(['L', 'M', 'H'], [[0], [0], [0]]))
        trs = dict(zip(['L', 'M', 'H'], [[0], [0], [0]]))          # key - condition: value - list of trials for each
        for t in trials:
            cond = CM['Target'][t]
            df = weights_to_edge_df(CM['A'][t])
            df = clean_df(df, 0.25)
            df.reset_index(drop=True, inplace=True)
            N = Nx.from_pandas_edgelist(df, 'Source', 'Target', 'Weight')  # network - use to get number of nodes
            if len(N.nodes) == 50:                                # if all nodes are present
                le = networkx.local_efficiency(N)                 # local efficiency
                ge = networkx.global_efficiency(N)                # global efficiency
                LE[cond].append(le)                               # add efficiencies to list according to condition
                GE[cond].append(ge)
                trs[cond].append(trs[cond][-1] + 1)               # running count of the number of trials for cond
        matrices = [LE, GE]
        measures = ['Local Efficiency', 'Global Efficiency']
        fig = plt.figure(figsize=[9, 7.5])                          # make a figure
        plt.xlabel('Trial', labelpad=20)                          # common xlabel
        plt.ylabel('Efficiency', labelpad=35)                     # common ylabel
        plt.title(f'Shift in Local and Global Efficiency', pad=25)  # common title
        plt.xticks([])                                            # suppress common axes ticks
        plt.yticks([])
        for condition, cond in cond_map.items():                  # for each condition
            ax = fig.add_subplot(3, 1, cond)                      # add a subplot to the figure
            ax.title.set_text(f'Condition {condition}')           # subplot title
            # ax.set_xticks(np.linspace(1, trs[condition][-1], num=len(trs['M'][1:]),
            # endpoint=True, dtype=int))                         # subplot xticks are trials for that condition
            for mat, meas in zip(matrices, measures):             # for each efficiency measure
                plt.plot(trs[condition][1:], mat[condition][1:], color=cond_color[condition],
                         label=meas)    # plot the efficiency over the relevant trials
        axes = fig.get_axes()                                     # force y axes to have same limits
        axes[0].get_shared_y_axes().join(*axes)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.5, top=0.85)                 # space for titles
        plt.legend()
        plt.show()
    else:                                                         # without considering condition
        LE = list()                                               # list of each efficiency measure
        GE = list()
        trs = list()                                              # list of trials with all nodes
        for t in trials:
            df = weights_to_edge_df(CM['A'][t])
            df = clean_df(df, 0.25)
            df.reset_index(drop=True, inplace=True)
            N = Nx.from_pandas_edgelist(df, 'Source', 'Target', 'Weight')
            if len(N.nodes) == 50:                                # if all nodes are present given the threshold
                le = networkx.local_efficiency(N)
                ge = networkx.global_efficiency(N)
                LE.append(le)                                     # add the efficiencies to their list
                GE.append(ge)
                trs.append(t)                                     # add the trial to the list
        matrices = [LE, GE]
        measures = ['Local Efficiency', 'Global Efficiency']
        plt.figure()                                              # make a figure
        plt.title(f'Shift in Local and Global Efficiency Over {len(trs)} Trials')
        plt.xlabel('Trial')
        plt.ylabel(f'Efficiency')
        for mat, meas in zip(matrices, measures):                 # for each efficiency measure
            # plot the efficiency over the relevant trials and label by efficiency measure
            plt.plot(trs, mat, label=meas)
        plt.legend()
        plt.show()


def avg_clustering_coeff_comparison(trials, condition_split=False):
    """ Plots the average clustering coefficient over trials with w/ or w/o splitting by condition
        Inputs:
          trials (list) - the trials that should be considered
          condition_split (bool) - split the plots by condition or not
        """
    CM = read_connectivity.get_CM()
    cond_map = {'L': 1, 'M': 2, 'H': 3}                                     # condition selector
    cond_color = {'L': 'red', 'M': 'yellow', 'H': 'green'}                  # condition color map
    if condition_split:
        # key - condition: value - average clustering coefficient for any amount of trials
        # ACC = Average Clustering Coefficient
        ACC = dict(zip(['L', 'M', 'H'], [[0], [0], [0]]))
        trs = dict(zip(['L', 'M', 'H'], [[0], [0], [0]]))
        for t in trials:
            cond = CM['Target'][t]
            df = weights_to_edge_df(CM['A'][t])
            df = clean_df(df, 0.25)
            df.reset_index(drop=True, inplace=True)
            N = Nx.from_pandas_edgelist(df, 'Source', 'Target', 'Weight')  # network - use to set attributes
            if len(N.nodes) == 50:                                         # if all nodes are present
                acc = networkx.average_clustering(N, weight='Weight')
                ACC[cond].append(acc)                                 # add the acc to the list for the proper condition
                trs[cond].append(trs[cond][-1] + 1)                   # running count of the number of trials for cond
        plt.figure()                                                  # make a figure
        plt.xlabel('Trial')
        plt.ylabel('Average Clustering Coefficient')
        plt.title(f'Shift in Average Clustering Coefficient by Condition')
        for condition, cond in cond_map.items():                           # for each condition
            plt.plot(trs[condition][1:], ACC[condition][1:], color=cond_color[condition],
                     label=condition)   # make a line of acc over applicable number of trials and label by condition
        plt.legend()
        plt.show()
    else:                                                                  # w/o considering condition
        ACC = list()                                                       # list of acc for all trials with all nodes
        trs = list()                                                       # list of all trials with all nodes
        for t in trials:
            df = weights_to_edge_df(CM['A'][t])
            df = clean_df(df, 0.25)
            df.reset_index(drop=True, inplace=True)
            N = Nx.from_pandas_edgelist(df, 'Source', 'Target', 'Weight')
            if len(N.nodes) == 50:                                         # if all nodes are present for a given threshold
                acc = networkx.average_clustering(N, weight='Weight')
                ACC.append(acc)                                            # add the weighted acc to the list
                trs.append(t)
        plt.figure()                                                       # make a figure
        plt.title(f'Shift in Average Clustering Coefficient Over {len(trs)} Trials')
        plt.xlabel('Trial')
        plt.ylabel('Average Clustering Coefficient')
        plt.plot(trs, ACC)                                                # plot all acc's over trials with all nodes
        plt.show()


def small_world_coeff_comparison(trials, condition_split=False):
    """ Plots the sigma small world coefficient over trials with w/ or w/o splitting by condition
            Inputs:
              trials (list) - the trials that should be considered
              condition_split (bool) - split the plots by condition or not
            """
    CM = read_connectivity.get_CM()
    cond_map = {'L': 1, 'M': 2, 'H': 3}                             # condition selector
    cond_color = {'L': 'red', 'M': 'yellow', 'H': 'green'}          # condition color map
    if condition_split:
        # key - condition: value - small world coefficient for any amount of trials
        # SWC = small world coefficient
        SWC = dict(zip(['L', 'M', 'H'], [[0], [0], [0]]))
        trs = dict(zip(['L', 'M', 'H'], [[0], [0], [0]]))
        for t in trials:
            cond = CM['Target'][t]
            df = weights_to_edge_df(CM['A'][t])
            df = clean_df(df, 0.25)
            df.reset_index(drop=True, inplace=True)
            N = Nx.from_pandas_edgelist(df, 'Source', 'Target', 'Weight')  # network - use to set attributes
            if len(N.nodes) == 50:                                  # if all nodes are present
                swc = networkx.sigma(N, niter=10)
                SWC[cond].append(swc)                               # add the swc to the list for the proper condition
                trs[cond].append(trs[cond][-1] + 1)                 # running count of the number of trials for cond
        plt.figure()                                                # make a figure
        plt.xlabel('Trial')
        plt.ylabel('Sigma')
        plt.title(f'Shift in Sigma Small World Coefficient by Condition')
        for condition, cond in cond_map.items():                    # for each condition
            plt.plot(trs[condition][1:], SWC[condition][1:], color=cond_color[condition],
                     label=condition)  # make a line of swc over applicable number of trials and label by condition
        plt.legend()
        print(trs)
        print(SWC)
        plt.axhline(1, linestyle='--', color='k')
        plt.show()
    else:                                                          # w/o considering condition
        SWC = list()                                               # list of swc for all trials with all nodes
        trs = list()                                               # list of all trials with all nodes
        for t in trials:
            df = weights_to_edge_df(CM['A'][t])
            df = clean_df(df, 0.25)
            df.reset_index(drop=True, inplace=True)
            N = Nx.from_pandas_edgelist(df, 'Source', 'Target', 'Weight')
            if len(N.nodes) == 50:                                 # if all nodes are present for a given threshold
                swc = networkx.sigma(N)
                SWC.append(swc)                                    # add the swc to the list
                trs.append(t)
        print(trs)
        print(SWC)
        plt.figure()                                               # make a figure
        plt.title(f'Shift in Sigma Small World Coefficient Over {len(trs)} Trials')
        plt.axhline(1, linestyle='--', color='k')
        plt.xlabel('Trial')
        plt.ylabel('Sigma')
        plt.plot(trs, SWC)                                         # plot all swc's over trials with all nodes
        plt.show()


def num_edges(trials, threshold=0.25):
    CM = read_connectivity.get_CM()
    edges = []
    self_edges = []
    trs = []
    for t in trials:
        df = weights_to_edge_df(CM['A'][t])
        df = clean_df(df, threshold)
        dfself = df[(df['Source'] - df['Target'] == 0)]  # subset of edges that are self-affecting
        df = df[(df['Source'] - df['Target'] != 0)]  # subset of edges that are not self-affecting
        df.reset_index(drop=True, inplace=True)
        dfself.reset_index(drop=True, inplace=True)
        edges.append(len(df))
        self_edges.append(len(dfself))
        trs.append(t)
    plt.figure()
    plt.plot(trs, edges, label='Two-Node Edges')
    plt.plot(trs, self_edges, label='Self-Affecting Edges')
    plt.axhline(50, linestyle='--', color='k')                      # maximum number of self-affecting edges
    plt.legend()
    plt.title(f'Number of Edges with Weight > {threshold} over {len(trs)} trials')
    plt.xlabel('Trial')
    plt.ylabel('Number of Edges')
    plt.show()

# For Debugging

# nodes = [1, 2, 3, 4, 5]
# A = np.array([[1, 4, 2, -5, 1],[-2, 4, 2, 5, 1],[3, -4, 2, 5, 1],[4, 4, -2, 5, 1],[5, 4, -2, 5, 1]])
# df = weights_to_edge_df(A)
# df = clean_df(df, 2)
CM = read_connectivity.get_CM()
tr = 3
df = weights_to_edge_df(CM['A'][tr])
df = clean_df(df, 0.2)
df.reset_index(drop=True, inplace=True)

N = Nx.from_pandas_edgelist(df, 'Source', 'Target', 'Weight', create_using=Nx.DiGraph()) # network - use to set attributes
network_graph = from_networkx(N, Nx.circular_layout, scale=20, center=(0, 0))  # graph - use to change visualization
plot = figure(tools="pan,wheel_zoom,save,reset", active_scroll='wheel_zoom', x_range=Range1d(-20.1, 20.1),
              y_range=Range1d(-20.1, 20.1), title='Test Network')
plot = arrow_renderer(df, N, network_graph, plot)
#in_degrees = signed_in_degree(df)
#out_degrees = signed_out_degree(df)
#connectivity_matrix(CM['A'][tr], tr, CM['Target'][tr], node_regions)
#circlelay = get_circle_layout(N, 20, (0, 0))
trs = list(range(22, 41))
# degree_comparison(trs, condition_split=True, region_group=True)
#efficiency_comparison(trs, False)
#avg_clustering_coeff_comparison(trs, False)
small_world_coeff_comparison(trs, True)
#num_edges(trs, 0.30)
