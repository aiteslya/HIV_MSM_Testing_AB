function distr=create_cas_dur_distr(varargin)

% this function creates numeric realization of the distribution of the
% duration of casual partnerships depending on what was passed in varargin
if nargin==1 % produce geometric distribution
    year=365;
    p=varargin{1}; % probability
    max_days=ceil(year/12);
    distr=zeros(1,max_days);
    for k=1:1:max_days
        distr(1,k)=p*(1-p)^(k-1);
    end
else
    myname=mfilename; 
    error([myname,': only geometric distribution has been realized']);
end

end