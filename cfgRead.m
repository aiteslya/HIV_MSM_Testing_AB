function cfg = cfgRead(filename) % read configuration file name and decode the configuration structure
	jscfg = fileread(filename);

	cfg = jsondecode(jscfg);
end
