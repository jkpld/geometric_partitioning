function A = accumLines(x1,x2,img_size,conn)
% ACCUMLINES ray trace a set of lines and accumulate all of the pixels for
% each line.
%
% x1 : N x 2 array, start points for a line
% x2 : N x 2 array, end points for a line
% img_size : the size of the image to which we are accumulating the points.
%
% x1(i,1) : coordinate of the 1st (!) dimension for the i'th line
% x1(i,2) : coordinate of the 2nd (!) dimension for the i'th line
%
% Note that x1(i,:) is NOT [x,y], but [y,x] as y is the first dimension of
% an image.

% Trace a line through a 2d grid and find each pixel the line intersects
% http://playtechs.blogspot.com/2007/03/raytracing-on-grid.html
% http://stackoverflow.com/questions/32328179/opencv-3-0-python-lineiterator
N = size(x1,1);
A = zeros(img_size);

if conn == 4
    
    for pt = 1:N
        
        dy = (x1(pt,1)-x2(pt,1));
        dx = (x1(pt,2)-x2(pt,2));
        y = [x1(pt,1), x2(pt,1)];
        x = [x1(pt,2), x2(pt,2)];
        
        if dy == 0 % horizontal line
            if x(2) > x(1)
                xs = x(1):x(2);
            else
                xs = x(1):-1:x(2);
            end
            ys = y(1)*ones(1,numel(xs));
            w = zeros(1,numel(xs));
        elseif dx == 0 % vertical line
            if y(2) > y(1)
                ys = y(1):y(2);
            else
                ys = y(1):-1:y(2);
            end
            xs = x(1)*ones(1,numel(ys));
            w = ones(1,numel(ys));
        else % sloped line
            if abs(dy) > abs(dx)
                slope = dx/dy;
                if y(2) > y(1)
                    ys = y(1):y(2);
                else
                    ys = y(1):-1:y(2);
                end
                xs = round(slope*(ys - y(1))) + x(1);
                w = ones(1,numel(ys));
            else
                slope = dy/dx;
                if x(2) > x(1)
                    xs = x(1):x(2);
                else
                    xs = x(1):-1:x(2);
                end
                ys = round(slope*(xs - x(1))) + y(1);
                w = abs(slope)*ones(1,numel(xs));
            end
        end
        
        inds = ys + (xs - 1)*img_size(1);
        A(inds) = A(inds) + w;
    end
    
elseif conn == 8
    
    for pt = 1:N
        
        dy = abs(x1(pt,1)-x2(pt,1));
        dx = abs(x1(pt,2)-x2(pt,2));
        y = x1(pt,1);
        x = x1(pt,2);
        n = 1 + dy + dx;
        if x1(pt,1) > x2(pt,1)
            y_inc = -1;
        else
            y_inc = 1;
        end
        
        if x1(pt,2) > x2(pt,2)
            x_inc = -1;
        else
            x_inc = 1;
        end
        
        err = dy - dx;
        
        dy = 2*dy; % Multipling by 2 makes everything in this problem an integer.
        dx = 2*dx;
        
        inds = zeros(n,1); % linear index
        
        i = 1;
        while i < n+1
            inds(i) = y + (x-1)*img_size(1);
            
            if err > 0
                y = y + y_inc;
                err = err - dx;
            else
                x = x + x_inc;
                err = err + dy;
            end
            i = i + 1;
        end
        
        A(inds) = A(inds) + 1;
    end
    
else
    error('accumLines:badconnectivity','connectivity must be 4 or 8')
end

end
