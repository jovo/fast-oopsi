function [lines splitmask] = CTgetLineSegments(im)
[m n] = size(im);
conn = CTgetConnections(im);

splitmask = zeros(m,n,'uint8');
pickedmask = zeros(m,n,'uint8');
splitmask(conn > 2) = 1;

%first desect lines
lines = {};
nLines = 0;
[x y] = find(conn == 1);
while length(x) > 0
    while length(x) >0
        x = x(1); y= y(1);
        nLines = nLines + 1;

        lines{nLines}.points = [x y];
        pickedmask(x,y) = 1;

        notdoneline = 1;
        while notdoneline > 0
            conn(x,y) = 0;
            im(x,y) = 0;

            %now get the next point
            temp = conn(x-1:x+1,y-1:y+1);
            [deltax deltay] = find(temp > 0); 
            if ~isempty(deltax)
                nextx = x + deltax(1) - 2;
                nexty = y + deltay(1) - 2;

                %end point
                if splitmask(nextx, nexty) > 0 
                    notdoneline = 0;
                else
                   if conn(nextx, nexty) == 2 %intermediate point
                       x = nextx; y = nexty; 
                   else
                       conn(nextx, nexty) = 0;
                       im(nextx,nexty) = 0;
                       notdoneline = 0;
                   end 
                end
                lines{nLines}.points = [lines{nLines}.points; [nextx nexty]];
                pickedmask(nextx,nexty) = 1;
            else
                notdoneline = 0;
            end
        end
        conn = CTgetConnections(im);
        [x y] = find(conn == 1);
    end
    im(pickedmask > 0) = 0;
    conn = CTgetConnections(im);
    [x y] = find(conn == 1);
end

%now deal with single points
im(pickedmask > 0) = 0;
[x y] = find(im > 0);
for i = 1:length(x)
    nLines = nLines + 1;
    lines{nLines}.points = [x(i) y(i)];
end
