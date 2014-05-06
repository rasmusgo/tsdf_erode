function [ img ] = draw_polygon( img, P, color )

% Draw polygon segments as lines
for i = 1:(size(P,2)-1)
    img = draw_line( img, P(:,i), P(:,i+1), color );
end
% Close the loop from end to start
img = draw_line( img, P(:,end), P(:,1), color );

end

