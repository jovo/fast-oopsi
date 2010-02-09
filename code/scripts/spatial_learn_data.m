iPath = '/Users/echidna/Desktop/roi.tif'; % path to tif stack
Nc    = 2;                                % number of neurons visible
V.dt  = 1/33.1;                           % time step size

% load up the image
if 0
    info    = imfinfo(iPath);
    nFrames = length(info);
    width   = info(1).Width;
    height  = info(1).Height;
else
    nFrames = 3816;
    width   = 21;
    height  = 21;
end

F = zeros(width*height,nFrames);

for ii = 1:nFrames
    F(:,ii) = reshape(imread(iPath,ii),width*height,1);
end

% set variables 
V.T      = nFrames;                     % # of time steps
V.Npixels= width*height;                % # of pixels in each image
V.Ncells = Nc;                          % # cells per ROI
V.h      = height;                      % height of frame (pixels)
V.w      = width;                       % width of frame (pixels)

% initialize parameters for plotting results after each pseudo-EM step
V.fast_thresh  = 1;                     % whether to threshold spike train before estimating 'a' and 'b' (we always keep this on)
V.fast_plot    = 1;                     % whether to plot filter with each iteration
V.fast_iter_max= 3;                     % # iterations of EM to estimate params
V.fast_ignore_post=1;

% estimate the spatial filters ftw
P.a = [mean(F(1:220,:));mean(F(221:end,:))]; %pca_approx(F',V.Ncells);
P.a = P.a';
P.b = 0.01 * mean(P.a);

% initialize other parameters
P.sig   = max(P.a(:))*1;               % stan dev of noise (indep for each pixel)
C_0     = 0;                           % initial calcium
tau     = rand(V.Ncells,1)/2+.05;        % decay time constant for each cell
P.gam   = 1-V.dt./tau(1:Nc);             % set gam
%P.lam   = P.sig;
lam     = linspace(0,4*pi,V.T);
P.lam = repmat(exp(1.6*(sin(lam))'),1,V.Ncells); % rate-ish, ie, lam*dt=# spikes per second
for j=2:V.Ncells
    P.lam(:,j) = exp(1.6*(sin(lam+pi))');
end
P.lam   = P.lam-repmat(min(P.lam),V.T,1);

fast_oopsi(F,V,P);