import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import matplotlib.image as image
from datetime import datetime
import sys
import pandas as pd
import os
  
  
class TrajectoryArtist:
    """Tool to manually draw or trace an image to make a 2D trajectory .csv file for motion phantom use"""
    def __init__(self):
        self.xdata = []
        self.ydata = []
        (self.figure, self.ax) = plt.subplots()
        # range of linear sliders
        self.ax.set_aspect(1)
        self.ax.set_xlim([0, 50])
        self.ax.set_ylim([0, 50])
        self.ax.set_xlabel("x (mm)")
        self.ax.set_ylabel("y (mm)")
        self.figure.canvas.mpl_connect('button_release_event', self.click_plot)
        self.figure.canvas.mpl_connect('key_release_event', self.keyboard_command)
  
    def start(self):
        """ Gets an imaage to trace if needed and displays the plot """
        self.get_background()
        print("Click to add a point, press r to remove last point, or press e to end trajectory and save points")
        plt.show()  # display the plot without a background
  
    def click_plot(self, event):
        """ Keeps track of the mouse click locations """
        self.draw_click(event)
        self.xdata.append(event.xdata)
        self.ydata.append(event.ydata)
  
    def draw_click(self, event):
        """ Draws the point of the click on the plot and adds it to the data list """
        oldptslines, = plt.plot(self.xdata+[event.xdata], self.ydata+[event.ydata], 'k')
        oldpts = plt.scatter(self.xdata, self.ydata, c='k')
        newpt = plt.scatter(x=event.xdata, y=event.ydata, c='r')
        # add to canvas
        event.canvas.figure.draw_artist(newpt)
        event.canvas.figure.gca().draw_artist(oldptslines)
        event.canvas.figure.gca().draw_artist(oldpts)
        event.canvas.draw()
        # print current data
        self.show_points(event)

    def get_background(self):
        """Get the file path of the background image"""
        self.image_path = input("Enter the full path of an image to trace or type n to draw free-hand: ")
        if self.image_path.upper() != "N" and os.path.exists(self.image_path):
            self.add_background()
        else:
            self.image_path = False

    def add_background(self):
        """ Adds background image to plot as a watermark """
        if not self.image_path:  # if there is no background image
            return
        with cbook.get_sample_data(self.image_path) as file:
            im = image.imread(file)
        # Add image as watermark to plot
        self.ax.imshow(im, extent=(10, 40, 10, 40), alpha=0.5)

    def keyboard_command(self, event):
        """ Called on a keyboard press event. Used to take user input as an interrupt """
        sys.stdout.flush()
        if event.key == 'r':  # remove last point
            x = self.xdata.pop()
            y = self.ydata.pop()
            self.reset_axes() # clear axes and put background image and labels back on
            # Replace points minus the last point
            oldptslines, = plt.plot(self.xdata, self.ydata, 'k')
            oldpts = plt.scatter(self.xdata, self.ydata, c='k')
            event.canvas.figure.gca().draw_artist(oldptslines)
            event.canvas.figure.gca().draw_artist(oldpts)
            event.canvas.draw()
            print("Removed Point: ({:.2f}, {:.2f})".format(x, y))
        elif event.key == 'e':  # exit artist mode
            plt.close(self.figure)
            self.save_points()
        sys.stdout.flush()

    def reset_axes(self):
        """ Removes all points and lines from axes """
        self.ax.clear()
        self.add_background()
        self.ax.set_aspect(1)
        self.ax.set_xlim([0, 50])
        self.ax.set_ylim([0, 50])
        self.ax.set_xlabel("x (mm)")
        self.ax.set_ylabel("y (mm)")

    def save_points(self):
        """ Saves points as a constant velocity-defined trajectory.csv file """
        self.save_path = input("Enter full path of where you want to save the trajectory file or press b to bypass saving: ")
        if self.save_path.upper() != "B" and os.path.exists(self.save_path):
            trajectory = {'x': self.xdata, 'y': self.ydata}
            df = pd.DataFrame.from_dict(trajectory)
            date_today = datetime.today().strftime("%m_%d_%y")
            current_time = datetime.now().strftime("%H_%M")
            full_path = os.path.join(self.save_path, "Trajectory_{}_{}.csv".format(date_today, current_time))
            df.to_csv(full_path, index=True)
            print("Saved as {}".format(full_path))
        else:
            print("Not saved")

    def show_points(self, event):
        """ Prints collected points """
        print("Trajectory Points (x, y)")
        for x, y in zip(self.xdata+[event.xdata], self.ydata+[event.ydata]):
            print("({:.2f}, {:.2f}) ".format(x, y))
  
ta = TrajectoryArtist()
ta.start()