
n = 1058;
imdata = imread('world.jpg');
s = imdata(:,:,2);
s = 1*(s<0);
s = reshape(s,[n,n]);
[u,sig,v] = svd(s);
num = 40;
imshow(u(:,1:num)*sig(1:num,1:num)*v(:,1:num)',[0 1])