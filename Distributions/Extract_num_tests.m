function add_tests = Extract_num_tests(sim_type, pars_num, traj_num)

    add_tests = zeros(pars_num, 1);

    for pars_counter = 1:pars_num
        add_tests(pars_counter, 1) = Extract_compar_interv_incid_pars(sim_type, pars_counter, traj_num);
    end

end