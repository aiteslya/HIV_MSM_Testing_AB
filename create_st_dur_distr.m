function distr=create_st_dur_distr(varargin)

% this function creates numeric realization of the distribution (pdf) of the
% duration of steady partnerships depending on what was passed in varargin


if nargin==1 % produce geometric distribution
    year=365;
    max_dur = 55; % maximum duration of a partnership in years

    p=varargin{1}; % probability
    max_days=max_dur*year;
    n = 1:max_days;
    distr = p * (1-p).^(n-1);
    distr = distr./sum(distr);
elseif nargin == 3 
    if isequal(varargin{1},'Erlang') | isequal(varargin{1},'erlang') % produce Erlang distribution
        year=365;
        max_dur = 55; % maximum duration of a partnership in years

        max_days = max_dur*year;
        mean_d = varargin{2};
        k = varargin{3};
        lambda = k/mean_d;
        x = 1:max_days;
        distr = (lambda^k).*(x.^(k-1)).*exp(-lambda*x)./(factorial(k-1));
        distr = distr./sum(distr);
    else
        myname=mfilename; 
        error([myname,': only geometric and Erlang distributions have been realized']);
    end
else
    myname=mfilename; 
    error([myname,': only geometric and Erlang distributions have been realized']);
end

end