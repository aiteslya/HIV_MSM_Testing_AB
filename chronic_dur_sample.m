function dur=chronic_dur_sample(varargin)
% this function samples a duration of a chronic stage of HIV infection
% without ART
% the details of the distribution to be sampled are contained in varargin

if nargin<1
    nameStr = mfilename;
    error([nameStr,': requires at least one input parameter']);
elseif nargin==1 % expect the mean of the distribution
    if ~isa(1,'double')
        nameStr = mfilename;
        error([nameStr,': requires at the input parameter to be double']);
    else
        % default sampling from exponential distribution
        lambda=1/varargin{1};
        dur=-log(1-rand)/lambda;
    end
elseif nargin==2
    mean_dur=varargin{1};
    lambda=1/mean_dur;
    distr_name=varargin{2};
    if isequal(distr_name,'exponential') | isequal(distr_name,'Exponential')
        dur=-log(1-rand)/lambda;
    else
        nameStr = mfilename;
        error([nameStr,': does not account for the ',distr_name,' distribution']);
    end
elseif nargin==3
    mean_dur=varargin{1};
    distr_name=varargin{2};
    shape_par=varargin{3};

    if isequal(distr_name,'erlang') | isequal(distr_name,'Erlang')
        scale_par=mean_dur/shape_par;
        dur=gamrnd(shape_par,scale_par);
    else
        nameStr = mfilename;
        error([nameStr,': does not account for the ',distr_name,' distribution']);
    end
else
    nameStr = mfilename;
    error([nameStr,': does not account for the number of parameters greater than 3']);
end
end