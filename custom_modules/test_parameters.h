/*
int number_of_permeating_solutes=3; //TODO: eventually fix oxygen issue (extra solute) so right now +1 permeating solutes since the 0th is oxygen
	int number_of_NP_solutes=1;
	std::string location="output/";
	int max_time_doub=300;//seconds
	double EG_pmv=0.0557414;
	double GLY_pmv;
	double TZP_k=0.05;//0.08;//0.045;//0.3;pbs//0.07;/0.09
	double max_TZP_length=10;//10 //Neighboorhood distance and breakage length
	double granulosa_k=0.3;//0.4;//0.4;
	double basal_k=0.6;//0.6;//0.6;
	double max_granulosa_connection_length=2; //Neighboorhood distance, gap junctions do not currently break
	double cell_compression_overlap=2;//how far membranes are allowed to overlap
	double Lp_oocyte=0.013;//0.0135;//0.0166666667;  //Combo is 0.0053//GLY is 0.007 //EG is 0.0166666666;//PBS is 0.0135;////um/atm/sec //EG 1/60=.016666667;https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4141568/;
	double Ps_oocyte=0.120;//0.12;//0.0; //.12;//0.0;um/sec//just refers to EG // combo is 0.798
	double Ps_oocyte_gly=0.107;//0.003;//0.107 //combo is 1.15
	double Lp_granulosa=0.00288776;//pbs=0.0114035; um/atm/s //parameters.doubles("granulosa_Lp");;0.00288776;
	double Ps_granulosa=0.0345;//um/s //parameters.doubles("granulosa_Ps");;//just refers to EG
	double Ps_granulosa_gly=0.0345;
	double final_external_osmolarity=2.185;//2.503;//.643;//.187;//.255;//1.6//2.185
	double holding_media_osmolarity=.2534;
	double loading_cpa_osmolarity=final_external_osmolarity-holding_media_osmolarity;
	double unloading_cpa_osmolarity=holding_media_osmolarity-final_external_osmolarity;
	double initial_granulosa_radius=5.0;//from config eventually
	double initial_oocyte_radius=58;//again should be global
	double follicle_radius =initial_oocyte_radius+(initial_granulosa_radius*2+initial_granulosa_radius*2+initial_granulosa_radius*2)+initial_granulosa_radius; // 3-4 layers
	double Temperature=296.65;//296.65;//23.5 Celsius=296.65

	std::vector <double> total_uptake;
	std::vector <std::vector <double> > all_total_uptakes;
	std::vector<Cell*> colored_cell_list;
	std::vector <int> number_of_cells_in_voxel;
	Cell* outtermost_granulosa;
	std::string max_time_str=std::to_string(max_time_doub);
	std::string CPA_name="EG";//"_15_percent_EG_dt_0p04";//should use to_string with actual parameters
	std::string test_name=max_time_str+"_"+CPA_name+"_";
	////////////////////////////////////vector parameters/////////////////////////////////////////////////////////////////


	std::vector<Cell*> basal_lamina_list; // list of cells connected to the BL at tissue creation
	double average_basement_distance=0;
	double number_of_basal_cells=0;
	//!setup textfile output////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::string output_TZP_distance_=location+test_name+ "_TZP_distance_TZP_k_" + std::to_string(int(tzp_k*10)) + ".csv";
	std::string output_Granulosa=location+test_name+ "_Granulosa.csv";
	std::string output_ecGranulosa=location+test_name+ "_external_concentration_Granulosa.csv";//not currently output
	std::string output_Oocyte=location+test_name+ "_Oocyte.csv";
	std::string output_ecOocyte=location+test_name+ "_external_concentration_Oocyte.csv";
	std::string output_Cell_0=location+test_name+ "_colored_cell_0.csv";
	std::string output_Cell_1=location+test_name+ "_colored_cell_1.csv";
	std::string output_Cell_2=location+test_name+ "_colored_cell_2.csv";
	std::string output_Cell_3=location+test_name+ "_colored_cell_3.csv";
	std::string output_ecCell_0=location+test_name+ "_external_concentration_cell_0.csv";
	std::string output_ecCell_1=location+test_name+ "_external_concentration_cell_1.csv";
	std::string output_ecCell_2=location+test_name+ "_external_concentration_cell_2.csv";
	std::string output_ecCell_3=location+test_name+ "_external_concentration_cell_3.csv";

	////////////////////////////////!test variables////////////
	std::vector<double> test_vector;
	std::vector<double> hole={10.0,120.0,0.0};
	double test1=0.0;
	double test2=0.0;
	double test3=0.0;
	*/
