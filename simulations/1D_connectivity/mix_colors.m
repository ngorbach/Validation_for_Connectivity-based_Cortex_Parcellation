function color_new = mix_colors(dist,color)

color_new = zeros(size(dist,1),3);
for k = 1:size(dist,2)
    color_new = color_new + 0.6*bsxfun(@times,color(k,:),dist(:,k));
end