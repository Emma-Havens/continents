function MARK2ELEM = mark_comm(ney, asp, L, MARK)

% Assumes structured, homogeneous grid

dx       = 2*L(1)/(asp*ney);
dy       = 2*L(2)/(ney);

x_ind    = 1 + floor((MARK.NODES(1,:)+L(1))./dx); % which column the marker belongs to
x_ind(x_ind>(asp*ney)) = asp*ney;

y_ind    = 1 + floor((MARK.NODES(2,:)+L(2))./dy); % which row the marker belongs to
y_ind(y_ind>ney) = ney;

MARK2ELEM = x_ind+(y_ind-1)*ney;
