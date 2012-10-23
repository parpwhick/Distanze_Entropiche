R = 29;
LX = 20;
LY = 20;
T = 201;

positions = zeros(R,T);
pos_labels = 1:LX;
for run = 1:R
filename = sprintf('run_%d.bin',run);
conf = fread(fopen(filename),inf,'int8');
conf = reshape(conf,LY,LX,T);

for t=1:T
    bubble_distribution = abs(sum(conf(:,:,t)));
    center = sum(bubble_distribution .* pos_labels) / sum(bubble_distribution);
    positions(run,t) = center;
end

end