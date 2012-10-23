function positions = plot_displacement(conf)

[LY,LX,T] = size(conf);

positions = zeros(T,1);
pos_labels = 1:LX;
for t=1:T
    bubble_distribution = abs(sum(conf(:,:,t)));
    center = sum(bubble_distribution .* pos_labels) / sum(bubble_distribution);
    positions(t) = center;
end

figure
plot(positions)
    