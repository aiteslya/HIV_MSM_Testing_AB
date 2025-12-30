function [scenarioName, pct1, pct2, pct3] = parseSimType(sim_type)
    %PARSESIMTYPE Extracts scenario name and numeric values from sim_type string
    % Handles special cases like 'ahi_XX_YY' and 'late_inc_WW_XX_YY_ZZ'
    
    tokens = split(sim_type, '_');
    n = numel(tokens);
    toFraction = @(t) str2double(t) ./ (10.^strlength(t));
    pct1 = [];
    pct2 = [];
    pct3 = [];

    % Special case for 'late_inc_WW_XX_YY_ZZ_AA_BB'
    if n == 6 && strcmp(tokens{1}, 'immigr') && strcmp(tokens{2}, 'late')
        t1 = tokens{3}; t2 = tokens{4}; t3 = tokens{5}; t4 = tokens{6};

        v1 = str2double(t1);
        v2 = str2double(t2);
        v3 = str2double(t3);
        v4 = str2double(t4);

        if all(~isnan([v1, v2, v3, v4]))
            scenarioName = 'immigr_late';
            pct1 = str2double(t1) + toFraction(t2); % WW.XX
            pct2 = str2double(t3) + toFraction(t4); % YY.ZZ
            return;
        end
    end

    % Special case for 'immigr_parts_WW_XX_YY_ZZ_AA_BB'

    if n == 8 && strcmp(tokens{1}, 'immigr') && strcmp(tokens{2}, 'parts')
        t1 = tokens{3}; t2 = tokens{4}; t3 = tokens{5}; t4 = tokens{6};
        t5 = tokens{7};  t6 = tokens{8}; 

        v1 = str2double(t1);
        v2 = str2double(t2);
        v3 = str2double(t3);
        v4 = str2double(t4);
        v5 = str2double(t5);
        v6 = str2double(t6);
        if all(~isnan([v1, v2, v3, v4, v5, v6]))
            scenarioName = 'immigr_parts';
            pct1 = str2double(t1) + toFraction(t2); % WW.XX
            pct2 = str2double(t3) + toFraction(t4); % YY.ZZ
            pct3 = str2double(t5) + toFraction(t6); % AA.BB
            return;
        end
    end

    % Special case for 'immigr_WW_XX'
    if n == 3 && strcmp(tokens{1}, 'immigr') 
        t1 = tokens{2}; t2 = tokens{3}; 

        v1 = str2double(t1);
        v2 = str2double(t2);

        if all(~isnan([v1, v2]))
            scenarioName = 'immigr';
            pct1 = str2double(t1) + toFraction(t2); % WW.XX
            return;
        end
    end

  % Special case for 'late_inc_WW_XX_YY_ZZ_AA_BB'
    if n == 8 && strcmp(tokens{1}, 'late') && strcmp(tokens{2}, 'inc')
        t1 = tokens{3}; t2 = tokens{4}; t3 = tokens{5}; t4 = tokens{6};
        t5 = tokens{7}; t6 = tokens{8};

        v1 = str2double(t1);
        v2 = str2double(t2);
        v3 = str2double(t3);
        v4 = str2double(t4);
        v5 = str2double(t5);
        v6 = str2double(t6);

        if all(~isnan([v1, v2, v3, v4, v5, v6]))
            scenarioName = 'late_inc';
            pct1 = str2double(t1) + toFraction(t2); % WW.XX
            pct2 = str2double(t3) + toFraction(t4); % YY.ZZ
            pct3 = str2double(t5) + toFraction(t6); % AA.BB
            return;
        end
    end


    % Special case for 'late_inc_WW_XX_YY_ZZ'
    if n == 6 && strcmp(tokens{1}, 'late') && strcmp(tokens{2}, 'inc')
        t1 = tokens{3}; t2 = tokens{4}; t3 = tokens{5}; t4 = tokens{6};

        v1 = str2double(t1);
        v2 = str2double(t2);
        v3 = str2double(t3);
        v4 = str2double(t4);
        if all(~isnan([v1, v2, v3, v4]))
            scenarioName = 'late_inc';
            pct1 = str2double(t1) + toFraction(t2); % WW.XX
            pct2 = str2double(t3) + toFraction(t4); % YY.ZZ
            return;
        end
    end

     % Special case for 'late_inc_WW_XX_YY_ZZ'
    if n == 6 && strcmp(tokens{1}, 'parts') && strcmp(tokens{2}, 'inc')
        t1 = tokens{3}; t2 = tokens{4}; t3 = tokens{5}; t4 = tokens{6};

        v1 = str2double(t1);
        v2 = str2double(t2);
        v3 = str2double(t3);
        v4 = str2double(t4);
        if all(~isnan([v1, v2, v3, v4]))
            scenarioName = 'parts_inc';
            pct1 = str2double(t1) + toFraction(t2); % WW.XX
            pct2 = str2double(t3) + toFraction(t4); % YY.ZZ
            return;
        end
    end

    % Special case for 'ahi_XX_YY'
    if n == 3 && strcmp(tokens{1}, 'ahi')
        t1 = tokens{2};
        t2 = tokens{3};
        if all(~isnan([str2double(t1), str2double(t2)]))
            scenarioName = 'ahi';
            pct1 = str2double(t1) + toFraction(t2);
            return;
        end
    end

    % General case with numeric suffixes
    if n >= 3
        t1 = tokens{end-1}; t2 = tokens{end};
        v1 = str2double(t1); v2 = str2double(t2);
        if ~isnan(v1) && ~isnan(v2)
            scenarioName = strjoin(tokens(1:end-2), '_');
            pct1 = toFraction(t1);
            pct2 = toFraction(t2);
            return;
        end
        v = str2double(t2);
        if ~isnan(v)
            scenarioName = strjoin(tokens(1:end-1), '_');
            pct1 = toFraction(t2);
            return;
        end
    end

    % Case with single numeric suffix
    if n == 2
        v = str2double(tokens{2});
        if ~isnan(v)
            scenarioName = tokens{1};
            pct1 = toFraction(tokens{2});
            return;
        end
    end

    % Fallback
    scenarioName = sim_type;

end
