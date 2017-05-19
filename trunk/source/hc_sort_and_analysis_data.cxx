// hc_sort_raw_data.cxx
// Standard libraries :
#include <iostream>
#include <bitset>
#include <fstream>

// - Bayeux/datatools:
#include <datatools/utils.h>
#include <datatools/io_factory.h>
#include <datatools/clhep_units.h>
// - Bayeux/geomtools:
#include <bayeux/geomtools/manager.h>
#include <bayeux/geomtools/id_mgr.h>
#include <bayeux/geomtools/id_selector.h>
// - Bayeux/mctools:
#include <mctools/simulated_data.h>
// - Bayeux/dpp:
#include <dpp/input_module.h>
#include <dpp/output_module.h>

// Falaise:
#include <falaise/falaise.h>
#include <falaise/snemo/geometry/gg_locator.h>
#include <falaise/snemo/geometry/calo_locator.h>

// Third part :
// Root :
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
// Boost :
#include "boost/filesystem.hpp"

int column_to_hc_half_zone(const int & column);

struct calo_hit_summary
{
  geomtools::geom_id geom_id;
  double energy = 0;
  double time = 0;
  geomtools::vector_3d left_most_hit_position; // The left most position
  bool geiger_association = false;
};

int main( int  argc_ , char **argv_  )
{
  falaise::initialize(argc_, argv_);
  int error_code = EXIT_SUCCESS;
  datatools::logger::priority logging = datatools::logger::PRIO_FATAL;

  // Parsing arguments
  int iarg = 1;
  bool is_input_file   = false;
  bool is_event_number = false;
  bool is_output_path  = false;
  bool is_display      = false;
  //  bool is_zone         = false;
  bool is_help         = false;

  std::string input_filename;
  std::string output_path;
  int arg_event_number  = -1;
  double calo_threshold = 50;
  //  int arg_hc_half_zone  = -1;

  while (iarg < argc_) {
    std::string arg = argv_[iarg];
    if (arg == "-i" || arg == "--input")
      {
	is_input_file = true;
	input_filename = argv_[++iarg];
      }

    else if (arg == "-op" || arg == "--output-path")
      {
	is_output_path = true;
	output_path = argv_[++iarg];
      }

    else if (arg == "-n" || arg == "--number")
      {
        is_event_number = true;
	arg_event_number    = atoi(argv_[++iarg]);
      }

    // else if (arg == "-z" || arg == "--zone")
    //   {
    //     is_zone = true;
    // 	arg_hc_half_zone  = atoi(argv_[++iarg]);
    //   }

    else if (arg == "-ct" || arg == "--calo-threshold")
      {
	calo_threshold = atoi(argv_[++iarg]);
      }

    else if (arg == "-d" || arg == "--display")
      {
	is_display = true;
      }

    else if (arg =="-h" || arg == "--help")
      {
	is_help = true;
      }

    iarg++;
  }

  if (is_help)
    {
      std::cerr << std::endl << "Goal : Sort raw data file to reduce the number of event / file " << std::endl << std::endl
		<< "Usage :" << std::endl << std::endl
		<< "$ BuildProducts/bin/half_commissioning-hc_sort_raw_data [OPTIONS] [ARGUMENTS]" << std::endl << std::endl
		<< "Allowed options: " << std::endl
		<< "-h  [ --help ]           produce help message" << std::endl
		<< "-i  [ --input ]          set an input file" << std::endl
		<< "-op [ --output path ]    set a path where all files are stored" << std::endl
		<< "-d  [ --display ]        display things for debug" << std::endl
		<< "-n  [ --number ]         set the number of events" << std::endl
		<< "Example : " << std::endl << std::endl;
      return 0;
    }

  try {
    std::clog << "INFO : Program sorting raw data file to reduce the number of event / file and produce a TTree in a ROOT file for analysis" << std::endl;

    std::string manager_config_file;
    manager_config_file = "@falaise:config/snemo/demonstrator/geometry/4.0/manager.conf";
    datatools::fetch_path_with_env(manager_config_file);
    datatools::properties manager_config;
    datatools::properties::read_config (manager_config_file,
					manager_config);
    geomtools::manager my_geom_manager;
    manager_config.update ("build_mapping", true);
    if (manager_config.has_key ("mapping.excluded_categories"))
      {
	manager_config.erase ("mapping.excluded_categories");
      }
    my_geom_manager.initialize (manager_config);

    std::string pipeline_simulated_data_filename = "";

    // Simulated Data "SD" bank label :
    std::string SD_bank_label = "SD";
    datatools::fetch_path_with_env(input_filename);
    if(is_input_file){
      pipeline_simulated_data_filename = input_filename;
    }else{
      // pipeline_simulated_data_filename = "/home/goliviero/software/Falaise/Users/GOliviero/HCAnalysis/trunk/data/Co60-100_row_0_column_0_SD.brio";
      pipeline_simulated_data_filename = "/sps/nemo/scratch/golivier/software/Falaise/Users/GOliviero/HCAnalysis/trunk/data/Co60-100_row_0_column_0_SD.brio";
    }
    std::clog << "INFO : INPUT FILENAME = " << input_filename << std::endl;
    datatools::fetch_path_with_env(pipeline_simulated_data_filename);


    // Vertex in the filename :
    std::string no_path_no_brio_filename = pipeline_simulated_data_filename;
    std::size_t found = 0;
    found = no_path_no_brio_filename.find_last_of("/");
    no_path_no_brio_filename.erase(0, found+1);
    found = no_path_no_brio_filename.find(".brio");
    no_path_no_brio_filename.erase(found, no_path_no_brio_filename.size());
    std::clog << "INFO : Filename without path and brio extension = " << no_path_no_brio_filename << std::endl;

    // Locators :
    int32_t my_module_number = 0;
    snemo::geometry::calo_locator calo_locator;
    calo_locator.set_geo_manager(my_geom_manager);
    calo_locator.set_module_number(my_module_number);
    calo_locator.initialize ();

    snemo::geometry::gg_locator gg_locator;
    gg_locator.set_geo_manager(my_geom_manager);
    gg_locator.set_module_number(my_module_number);
    gg_locator.initialize ();

    // Number of events :
    int event_number = -1;
    if (is_event_number) event_number = arg_event_number;
    else                 event_number = 10;

    std::clog << "INFO : Event number = " << event_number << std::endl;

    // Event reader :
    dpp::input_module reader;
    datatools::properties reader_config;
    reader_config.store ("logging.priority", "debug");
    reader_config.store ("max_record_total", event_number);
    reader_config.store ("files.mode", "single");
    reader_config.store ("files.single.filename", pipeline_simulated_data_filename);
    reader.initialize_standalone (reader_config);
    datatools::multi_properties iMetadataStore = reader.get_metadata_store();
    reader.tree_dump(std::clog, "Simulated data reader module");

    // Output path :
    if (is_output_path) datatools::fetch_path_with_env(output_path);
    //else output_path = "/home/goliviero/software/Falaise/Users/GOliviero/HCAnalysis/trunk/output_default/";
    // CC output path :
    else output_path = "/sps/nemo/scratch/golivier/software/Falaise/Users/GOliviero/HCAnalysis/trunk/output_default/";

    // Create output dir
    boost::filesystem::path rootPath(output_path);
    boost::system::error_code returnedError;
    boost::filesystem::create_directories(rootPath, returnedError);
    DT_THROW_IF(returnedError, std::logic_error, "Path : " + output_path + " was not created ! ");
    std::clog << "INFO : OUTPUT PATH = " << output_path << std::endl;

    std::string analysis_path = output_path + "analysis/";
    std::clog << "INFO : ANALYSIS PATH = " << analysis_path << std::endl;

    std::string sorted_path = output_path + "sorted_raw_data/";
    std::clog << "INFO : SORTED PATH = " << sorted_path << std::endl;
    
    std::string other_brio_files_path = analysis_path + "other_brio_files/";
    std::clog << "INFO : OTHER BRIO FILES PATH = " << other_brio_files_path << std::endl;
  
    boost::filesystem::path rootPath1(analysis_path);
    boost::filesystem::create_directories(rootPath1, returnedError);
    DT_THROW_IF(returnedError, std::logic_error, "Path : " + analysis_path + " was not created ! ");

    boost::filesystem::path rootPath2(sorted_path);
    boost::filesystem::create_directories(rootPath2, returnedError);
    DT_THROW_IF(returnedError, std::logic_error, "Path : " + sorted_path + " was not created ! ");

    boost::filesystem::path rootPath3(other_brio_files_path);
    boost::filesystem::create_directories(rootPath3, returnedError);
    DT_THROW_IF(returnedError, std::logic_error, "Path : " + other_brio_files_path + " was not created ! ");

    // Event record :
    datatools::things ER;

    /********************************************************************/
    /* DESIGN OF THE SOFT : JUST FOR ZONE 2 FOR THE MOMENT (no need to  */
    /* simulate all half zones (unefficient) maybe later. TO DO : adapt */
    /* the software for all half zones (take into account row/column or */
    /* free spots vertexes (row / column approach needs sometimes to be */
    /* process 2 times (column 3 4 5 count for zone 0 and 1 for example */
    /********************************************************************/

    unsigned int hc_half_zone_number = 2;

    // Selection rules for zone number = 2

    /* +-----+
       |     |......... GG cells #row 1
       |calo |......... GG cells #row 0
       |     |
       +-----+                          */

    std::size_t main_calo_column = hc_half_zone_number;
    std::size_t geiger_hc_zone_inf_limit = hc_half_zone_number * 6 - 3;
    std::size_t geiger_hc_zone_sup_limit = geiger_hc_zone_inf_limit + 5;

    // Id selector rules for half zone 2 :
    std::string hc_main_calo_half_zone_rules = "category='calorimeter_block' module={0} side={1} column={" + std::to_string(main_calo_column) + "} row={*} part={*}";
    // Only column 11 and 12 are sensitive because electronics for commissioning
    std::string hc_geiger_half_zone_rules = "category='drift_cell_core' module={0} side={1} layer={*} row={"
      + std::to_string(geiger_hc_zone_inf_limit + 2) + ";"
      + std::to_string(geiger_hc_zone_inf_limit + 3) + "}";

    geomtools::id_selector my_hc_main_calo_id_selector(my_geom_manager.get_id_mgr());
    my_hc_main_calo_id_selector.initialize(hc_main_calo_half_zone_rules);
    if (is_display) my_hc_main_calo_id_selector.dump(std::clog, "Main calo ID selector: ");


    geomtools::id_selector my_hc_geiger_id_selector(my_geom_manager.get_id_mgr());
    my_hc_geiger_id_selector.initialize(hc_geiger_half_zone_rules);
    if (is_display) my_hc_geiger_id_selector.dump(std::clog, "Geiger ID selector: ");

    //==============================================//
    //          output files  and writers           //
    //==============================================//

    // Name of sorted (matching rules) SD output file :
    std::string sorted_sd_brio = sorted_path + no_path_no_brio_filename + "_match_rules.brio";

    // Event writer for sorted SD :
    dpp::output_module sorted_writer;
    datatools::properties sorted_writer_config;
    sorted_writer_config.store ("logging.priority", "fatal");
    sorted_writer_config.store ("files.mode", "single");
    sorted_writer_config.store ("files.single.filename", sorted_sd_brio);
    sorted_writer.grab_metadata_store() = iMetadataStore;
    sorted_writer.initialize_standalone(sorted_writer_config);

    // Name of full track (9 layers hit) SD output file :
    std::string full_track_sd_brio = other_brio_files_path + no_path_no_brio_filename + "_match_rules_full_track.brio";

    // Event writer for full track (9 layers hit) :
    dpp::output_module full_track_writer;
    datatools::properties full_track_writer_config;
    full_track_writer_config.store ("logging.priority", "fatal");
    full_track_writer_config.store ("files.mode", "single");
    full_track_writer_config.store ("files.single.filename", full_track_sd_brio);
    full_track_writer.grab_metadata_store() = iMetadataStore;
    full_track_writer.initialize_standalone(full_track_writer_config);

    // Name of full track (9 layers hit) SD output file :
    std::string two_calos_one_full_track_sd_brio = other_brio_files_path + no_path_no_brio_filename + "_match_rules_two_calos_one_full_track.brio";

    // Event writer for 2 calos one full track (9 layers hit) :
    dpp::output_module two_calos_one_full_track_writer;
    datatools::properties two_calos_one_full_track_writer_config;
    two_calos_one_full_track_writer_config.store ("logging.priority", "fatal");
    two_calos_one_full_track_writer_config.store ("files.mode", "single");
    two_calos_one_full_track_writer_config.store ("files.single.filename", two_calos_one_full_track_sd_brio);
    two_calos_one_full_track_writer.grab_metadata_store() = iMetadataStore;
    two_calos_one_full_track_writer.initialize_standalone(two_calos_one_full_track_writer_config);

    // Output ROOT file :
    std::string string_buffer = analysis_path + no_path_no_brio_filename + ".root";
    datatools::fetch_path_with_env(string_buffer);

    TFile* root_file = new TFile(string_buffer.c_str(), "RECREATE");
    // TTree* root_tree = new TTree("HC_Analysis", "Half Commissioning analysis");

    TH1F * calo_number_TH1F = new TH1F(string_buffer.c_str(),
				       Form("Number of calorimeter touched, zone %i;", hc_half_zone_number),
				       13, 0, 13);
    string_buffer = "geiger_number_TH1F";
    TH1F * geiger_number_TH1F = new TH1F(string_buffer.c_str(),
					 Form("Number of Geiger cells touched, zone %i;", hc_half_zone_number),
					 54, 0, 54);
    string_buffer = "calo_total_energy_TH1F";
    TH1F * calo_total_energy_TH1F = new TH1F(string_buffer.c_str(),
					     Form("Calo total energy spectrum, zone %i;", hc_half_zone_number),
					     300, 0, 3000);
    string_buffer = "calo_total_time_TH1F";
    TH1F * calo_total_time_TH1F = new TH1F(string_buffer.c_str(),
					   Form("Time start for a calo, zone %i;", hc_half_zone_number),
					   150, 0, 15);
    string_buffer = "calo_total_distribution_TH2F";
    TH2F * calo_total_distribution_TH2F = new TH2F(string_buffer.c_str(),
						   Form("Calo total distribution, zone %i;", hc_half_zone_number),
						   20, 0, 20,
						   14, 0, 14);
    string_buffer = "geiger_total_distribution_TH2F";
    TH2F * geiger_total_distribution_TH2F = new TH2F(string_buffer.c_str(),
						     Form("Geiger total distribution, zone %i;", hc_half_zone_number),
						     7, geiger_hc_zone_inf_limit, geiger_hc_zone_sup_limit+2,
						     10, 0, 10);
    // One calo :
    string_buffer = "one_calo_energy_TH1F";
    TH1F * one_calo_energy_TH1F = new TH1F(string_buffer.c_str(),
					   Form("One calo energy, zone %i;", hc_half_zone_number),
					   150, 0, 1500);
    string_buffer = "one_calo_distribution_TH2F";
    TH2F * one_calo_distribution_TH2F = new TH2F(string_buffer.c_str(),
						 Form("One calo distribution, zone %i;", hc_half_zone_number),
						 20, 0, 20,
						 14, 0, 14);
    // One calo no track :
    string_buffer = "one_calo_no_track_energy_TH1F";
    TH1F * one_calo_no_track_energy_TH1F = new TH1F(string_buffer.c_str(),
						    Form("One calo no track energy, zone %i;", hc_half_zone_number),
						    150, 0, 1500);
    string_buffer = "one_calo_no_track_distribution_TH2F";
    TH2F * one_calo_no_track_distribution_TH2F = new TH2F(string_buffer.c_str(),
							  Form("One calo no track, zone %i;", hc_half_zone_number),
							  20, 0, 20,
							  14, 0, 14);
    // One calo one track :
    string_buffer = "one_calo_one_track_energy_TH1F";
    TH1F * one_calo_one_track_energy_TH1F = new TH1F(string_buffer.c_str(),
						     Form("One calo one track energy, zone %i;", hc_half_zone_number),
						     150, 0, 1500);
    string_buffer = "one_calo_one_track_distribution_TH2F";
    TH2F * one_calo_one_track_distribution_TH2F = new TH2F(string_buffer.c_str(),
							   Form("One calo one track, zone %i;", hc_half_zone_number),
							   20, 0, 20,
							   14, 0, 14);
    string_buffer = "one_calo_one_track_geiger_distribution_TH2F";
    TH2F * one_calo_one_track_geiger_distribution_TH2F = new TH2F(string_buffer.c_str(),
								  Form("One calo one track geiger distribution, zone %i;", hc_half_zone_number),
								  7, geiger_hc_zone_inf_limit, geiger_hc_zone_sup_limit+2,
								  10, 0, 10);

    // Two calos :
    string_buffer = "two_calo_energy_min_TH1F";
    TH1F * two_calo_energy_min_TH1F = new TH1F(string_buffer.c_str(),
					       Form("Two calos energy min spectrum, zone %i;", hc_half_zone_number),
					       150, 0, 1500);
    string_buffer = "two_calo_energy_max_TH1F";
    TH1F * two_calo_energy_max_TH1F = new TH1F(string_buffer.c_str(),
					       Form("Two calos energy max spectrum, zone %i;", hc_half_zone_number),
					       150, 0, 1500);
    string_buffer = "two_calo_total_energy_TH1F";
    TH1F * two_calo_total_energy_TH1F = new TH1F(string_buffer.c_str(),
						 Form("Two calos total energy spectrum, zone %i;", hc_half_zone_number),
						 300, 0, 3000);
    string_buffer = "two_calo_delta_time_TH1F";
    TH1F * two_calo_delta_time_TH1F = new TH1F(string_buffer.c_str(),
					       Form("Delta time (ns) between 2 calos, zone %i;", hc_half_zone_number),
					       100, 0, 15);
    string_buffer = "two_calo_distribution_TH2F";
    TH2F * two_calo_distribution_TH2F = new TH2F(string_buffer.c_str(),
						 Form("Two calos distrib,, zone %i;", hc_half_zone_number),
						 20, 0, 20,
						 14, 0, 14);
    string_buffer = "two_calo_geiger_distribution_TH2F";
    TH2F * two_calo_geiger_distribution_TH2F = new TH2F(string_buffer.c_str(),
							Form("Two calos geiger distribution, zone %i;", hc_half_zone_number),
							7, geiger_hc_zone_inf_limit, geiger_hc_zone_sup_limit+2,
							10, 0, 10);
    // Two calos no track :
    string_buffer = "two_calo_no_track_energy_min_TH1F";
    TH1F * two_calo_no_track_energy_min_TH1F = new TH1F(string_buffer.c_str(),
							Form("Two calos no track energy min spectrum, zone %i;", hc_half_zone_number),
							150, 0, 1500);
    string_buffer = "two_calo_no_track_energy_max_TH1F";
    TH1F * two_calo_no_track_energy_max_TH1F = new TH1F(string_buffer.c_str(),
							Form("Two calos no track energy max spectrum, zone %i;", hc_half_zone_number),
							150, 0, 1500);
    string_buffer = "two_calo_no_track_total_energy_TH1F";
    TH1F * two_calo_no_track_total_energy_TH1F = new TH1F(string_buffer.c_str(),
							  Form("Two calos no track total energy spectrum, zone %i;", hc_half_zone_number),
							  300, 0, 3000);
    string_buffer = "two_calo_no_track_delta_time_TH1F";
    TH1F * two_calo_no_track_delta_time_TH1F = new TH1F(string_buffer.c_str(),
							Form("Delta time (ns) between 2 calo no track, zone %i;", hc_half_zone_number),
							100, 0, 15);
    string_buffer = "two_calo_no_track_distribution_TH2F";
    TH2F * two_calo_no_track_distribution_TH2F = new TH2F(string_buffer.c_str(),
							  Form("Two calos no track distrib,, zone %i;", hc_half_zone_number),
							  20, 0, 20,
							  14, 0, 14);
    // Two calos one track :
    string_buffer = "two_calo_one_track_energy_min_TH1F";
    TH1F * two_calo_one_track_energy_min_TH1F = new TH1F(string_buffer.c_str(),
							 Form("Two calos one track energy min spectrum, zone %i;", hc_half_zone_number),
							 150, 0, 1500);
    string_buffer = "two_calo_one_track_energy_max_TH1F";
    TH1F * two_calo_one_track_energy_max_TH1F = new TH1F(string_buffer.c_str(),
							 Form("Two calos one track energy max spectrum, zone %i;", hc_half_zone_number),
							 150, 0, 1500);
    string_buffer = "two_calo_one_track_total_energy_TH1F";
    TH1F * two_calo_one_track_total_energy_TH1F = new TH1F(string_buffer.c_str(),
							   Form("Two calos one track total energy spectrum, zone %i;", hc_half_zone_number),
							   300, 0, 3000);
    string_buffer = "two_calo_one_track_electron_energy_TH1F";
    TH1F * two_calo_one_track_electron_energy_TH1F = new TH1F(string_buffer.c_str(),
							      Form("Two calos one track electron energy spectrum, zone %i;", hc_half_zone_number),
							      150, 0, 1500);
    string_buffer = "two_calo_one_track_gamma_energy_TH1F";
    TH1F * two_calo_one_track_gamma_energy_TH1F = new TH1F(string_buffer.c_str(),
							   Form("Two calos one track gamma energy spectrum, zone %i;", hc_half_zone_number),
							   150, 0, 1500);
    string_buffer = "two_calo_one_track_delta_time_TH1F";
    TH1F * two_calo_one_track_delta_time_TH1F = new TH1F(string_buffer.c_str(),
							 Form("Delta time (ns) between 2 calo one track, zone %i;", hc_half_zone_number),
							 100, 0, 15);
    string_buffer = "two_calo_one_track_angle_TH1F";
    TH1F * two_calo_one_track_angle_TH1F = new TH1F(string_buffer.c_str(),
						    Form("Angle (in degrees) between the vertex and the calorimeters hit, zone %i;", hc_half_zone_number),
						    90, 0, 90);
    string_buffer = "two_calo_one_track_calo_interaction_distribution_TH2F";
    TH2F * two_calo_one_track_calo_interaction_distribution_TH2F = new TH2F(string_buffer.c_str(),
									    Form("Two calos one track calo interaction distribution,, zone %i;", hc_half_zone_number),
									    10, -50, 150,
									    90, -1500, 1500);
    string_buffer = "two_calo_one_track_distribution_TH2F";
    TH2F * two_calo_one_track_distribution_TH2F = new TH2F(string_buffer.c_str(),
							   Form("Two calos one track distrib,, zone %i;", hc_half_zone_number),
							   20, 0, 20,
							   14, 0, 14);
    string_buffer = "two_calo_one_track_geiger_distribution_TH2F";
    TH2F * two_calo_one_track_geiger_distribution_TH2F = new TH2F(string_buffer.c_str(),
								  Form("Two calos one track geiger distribution, zone %i;", hc_half_zone_number),
								  7, geiger_hc_zone_inf_limit, geiger_hc_zone_sup_limit+2,
								  10, 0, 10);
    string_buffer = "calo_distribution_full_track_TH2F";
    TH2F * calo_distribution_full_track_TH2F = new TH2F(string_buffer.c_str(),
							Form("Calo full track distribution, zone %i;", hc_half_zone_number),
							20, 0, 20,
							14, 0, 14);

    // Event counter :
    int event_id    = 0;

    while (!reader.is_terminated())
      {
	if(is_display) std::clog << "INFO : Event #" << event_id << std::endl;
	reader.process(ER);

	bool match_rules_event = false;
	bool full_track_event  = false;

	// A plain `mctools::simulated_data' object is stored here :
	if (ER.has(SD_bank_label) && ER.is_a<mctools::simulated_data>(SD_bank_label))
	  {
	    // Access to the "SD" bank with a stored `mctools::simulated_data' :
	    const mctools::simulated_data & SD = ER.get<mctools::simulated_data>(SD_bank_label);

	    // New SD bank flaged for Geiger cell already hit
	    mctools::simulated_data flaged_SD = SD;

	    // First loop on all hits to see if ot match rules and is full track or not, also tag GG cells hit several times :

	    std::map<geomtools::geom_id, calo_hit_summary> calo_hit_map;
	    // If main calo hits :
	    if (flaged_SD.has_step_hits("calo"))
	      {
		// const size_t number_of_main_calo_hits = flaged_SD.get_number_of_step_hits("calo");
		mctools::simulated_data::hit_handle_collection_type BSHC = flaged_SD.get_step_hits("calo");
		if (is_display) std::clog << "BSCH calo step hits # = " << BSHC.size() << std::endl;

		for (mctools::simulated_data::hit_handle_collection_type::const_iterator i = BSHC.begin();
		     i != BSHC.end();
		     i++)
		  {
		    const mctools::base_step_hit & BSH = i->get();
		    // extract the corresponding geom ID:
		    const geomtools::geom_id & main_calo_gid = BSH.get_geom_id();
		    if (my_hc_main_calo_id_selector.match(main_calo_gid))
		      {
			match_rules_event = true;
			bool calo_hit_is_in_map = calo_hit_map.find(main_calo_gid) != calo_hit_map.end();
			// int column = main_calo_gid.get(2);
			// int row = main_calo_gid.get(3);

			if (!calo_hit_is_in_map)
			  {
			    // Create a new calo hit and insert in the map
			    calo_hit_summary new_calo_hit;
			    new_calo_hit.geom_id = main_calo_gid;
			    new_calo_hit.energy = BSH.get_energy_deposit();
			    new_calo_hit.time = BSH.get_time_start();
			    new_calo_hit.left_most_hit_position = BSH.get_position_start();
			    calo_hit_map.insert(std::pair<geomtools::geom_id, calo_hit_summary>(main_calo_gid, new_calo_hit));
			  }
			else
			  {
			    // Update the existing calo hit (add energy and check the time to be the t_start)min)
			    calo_hit_map.find(main_calo_gid)->second.energy+= BSH.get_energy_deposit();
			    if (BSH.get_time_start() < calo_hit_map.find(main_calo_gid)->second.time) calo_hit_map.find(main_calo_gid)->second.time = BSH.get_time_start();

			    if (BSH.get_position_start().getX() < calo_hit_map.find(main_calo_gid)->second.left_most_hit_position.getX())
			      {
				calo_hit_map.find(main_calo_gid)->second.left_most_hit_position = BSH.get_position_start();
			      }
			  }
		      } // end of match rules
		  } // end of for i BSHC

		// Check if all calo summary hit pass the threshold :
		for (std::map<geomtools::geom_id, calo_hit_summary>::iterator it_calo = calo_hit_map.begin();
		     it_calo != calo_hit_map.end();)
		  {
		    if (it_calo->second.energy * 1000 < calo_threshold)
		      {
			it_calo = calo_hit_map.erase(it_calo);
		      }
		    else it_calo++;
		  }

	      } // end of if has step hits "calo"

	    std::set<geomtools::geom_id> geiger_hit_set;
	    std::vector<geomtools::vector_3d> collection_position_last_geiger_hit;
	    std::size_t geiger_last_layer = 8;

	    if (flaged_SD.has_step_hits("gg"))
	      {
		const size_t number_of_gg_hits = flaged_SD.get_number_of_step_hits("gg");

		// We have to flag the gg cells already hit before (maybe take into account the dead time of a GG cell)
		for (size_t ihit = 0; ihit < number_of_gg_hits; ihit++)
		  {
		    mctools::base_step_hit & geiger_hit = flaged_SD.grab_step_hit("gg", ihit);
		    for (size_t jhit = ihit + 1; jhit < number_of_gg_hits; jhit++)
		      {
			mctools::base_step_hit & other_geiger_hit = flaged_SD.grab_step_hit("gg", jhit);
			if (geiger_hit.get_geom_id() == other_geiger_hit.get_geom_id())
			  {
			    const double gg_hit_time       = geiger_hit.get_time_start();
			    const double other_gg_hit_time = other_geiger_hit.get_time_start();
			    if (gg_hit_time > other_gg_hit_time)
			      {
				bool geiger_already_hit = true;
				if (!geiger_hit.get_auxiliaries().has_flag("geiger_already_hit")) geiger_hit.grab_auxiliaries().store("geiger_already_hit", geiger_already_hit);
			      }
			    else
			      {
				bool geiger_already_hit = true;
				if (!other_geiger_hit.get_auxiliaries().has_flag("geiger_already_hit")) other_geiger_hit.grab_auxiliaries().store("geiger_already_hit", geiger_already_hit);
			      }
			  } // end of if get_geom_id
		      } // end of jhit
		  } // end of ihit

		mctools::simulated_data::hit_handle_collection_type BSHC_gg = flaged_SD.get_step_hits("gg");
		if (is_display) std::clog << "BSCH geiger step hits # = " << BSHC_gg.size() << std::endl;
		for (mctools::simulated_data::hit_handle_collection_type::const_iterator i = BSHC_gg.begin();
		     i != BSHC_gg.end();
		     i++)
		  {
		    const mctools::base_step_hit & BSH = i->get();
		    // if (is_display) BSH.tree_dump(std::clog, "A Geiger Base Step Hit : ", "INFO : ");

		    // Second loop on Geiger hits, ignore cells hit 2 or more times.
		    if (BSH.get_auxiliaries().has_flag("geiger_already_hit") || BSH.get_auxiliaries().has_flag("other_geiger_already_hit")) {}
		    else
		      {
			const geomtools::geom_id & geiger_gid = BSH.get_geom_id();
			if (my_hc_geiger_id_selector.match(geiger_gid))
			  {
			    geiger_hit_set.insert(geiger_gid);
			    match_rules_event = true;
			    // If the last Geiger at layer 8 is hit, push back position for calorimeter association
			    if (geiger_gid.get(2) == geiger_last_layer) // last layer
			      {
				collection_position_last_geiger_hit.push_back(BSH.get_position_stop());
			      }
			  }
		      }
		  } // end of for
	      } // end of if has step hits "gg"

	    /*******************************************************/
	    /* Begining analysis (based on calo map and geiger set */
	    /*******************************************************/

	    // Fill histograms in ROOT file :
	    if (match_rules_event)
	      {
		sorted_writer.process(ER);

		int number_of_main_calo = calo_hit_map.size();
		int number_of_geiger = geiger_hit_set.size();

		calo_number_TH1F->Fill(number_of_main_calo);
		geiger_number_TH1F->Fill(number_of_geiger);

		// For each calorimeter, add it in the histogram
		for (std::map<geomtools::geom_id, calo_hit_summary>::const_iterator it_calo = calo_hit_map.begin(); it_calo != calo_hit_map.end(); it_calo++)
		  {
		    calo_total_energy_TH1F->Fill(it_calo->second.energy * 1000);
		    calo_total_time_TH1F->Fill(it_calo->second.time);
		    int column = it_calo->first.get(2);
		    int row = it_calo->first.get(3);
		    calo_total_distribution_TH2F->Fill(column, row);
		  }

		std::bitset<9> layer_projection = 0x0;

		// For each Geiger cell, add it in the histogram
		for (std::set<geomtools::geom_id>::const_iterator it_geiger = geiger_hit_set.begin(); it_geiger != geiger_hit_set.end(); it_geiger++)
		  {
		    int layer = it_geiger->get(2);
		    int row   = it_geiger->get(3);
		    // if (is_display) std::clog << "Layer  = " << layer << " Row = " << row << std::endl;
		    geiger_total_distribution_TH2F->Fill(row, layer);
		    layer_projection.set(layer, true);
		  }
		int number_of_layer = layer_projection.count();
		if (number_of_layer == 9) full_track_event = true;

		if (number_of_main_calo > 0)
		  {
		    // If 1+ calo & Full track
		    if (full_track_event)
		      {
			for (std::map<geomtools::geom_id, calo_hit_summary>::const_iterator it_calo = calo_hit_map.begin(); it_calo != calo_hit_map.end(); it_calo++)
			  {
			    int column = it_calo->first.get(2);
			    int row = it_calo->first.get(3);
			    calo_distribution_full_track_TH2F->Fill(column, row);
			  }
			full_track_writer.process(ER);
		      }

		    // One calo only :
		    if (number_of_main_calo == 1)
		      {
			one_calo_energy_TH1F->Fill(calo_hit_map.begin()->second.energy * 1000);
			for (std::map<geomtools::geom_id, calo_hit_summary>::const_iterator it_calo = calo_hit_map.begin(); it_calo != calo_hit_map.end(); it_calo++)
			  {
			    int column = it_calo->first.get(2);
			    int row = it_calo->first.get(3);
			    one_calo_distribution_TH2F->Fill(column, row);
			  }

			// One calo, no track at all
			if (number_of_layer == 0)
			  {
			    one_calo_no_track_energy_TH1F->Fill(calo_hit_map.begin()->second.energy * 1000);
			    for (std::map<geomtools::geom_id, calo_hit_summary>::const_iterator it_calo = calo_hit_map.begin(); it_calo != calo_hit_map.end(); it_calo++)
			      {
				int column = it_calo->first.get(2);
				int row = it_calo->first.get(3);
				one_calo_no_track_distribution_TH2F->Fill(column, row);
			      }
			  }

			// One calo, one full track (9 layers)
			if (number_of_layer == 9)
			  {
			    // Check the Z of calo and last layer geiger
			    calo_hit_summary calo_1 = calo_hit_map.begin()->second;
			    int side = 1;
			    // int column_1 = calo_1.geom_id.get(2);
			    int row_1 = calo_1.geom_id.get(3);

			    double calo_1_zposition = calo_locator.get_row_z(side, row_1);
			    if (is_display) std::clog << "Calo Z pos : " << calo_1_zposition << std::endl;
			    for (unsigned int i = 0; i < collection_position_last_geiger_hit.size(); i++)
			      {
				if (is_display) std::clog << "GG Z pos : " << collection_position_last_geiger_hit[i].z() << std::endl;
				if (calo_1_zposition > collection_position_last_geiger_hit[i].z() - 185
				    && calo_1_zposition < collection_position_last_geiger_hit[i].z() + 185)
				  {
				    if (is_display)
				      {
					std::clog << "Calorimeter and Geiger association" << std::endl;
					std::clog << "Calo 1 Z: " << calo_1_zposition << std::endl;
					std::clog << "GG hit #"<< i << ' ' << collection_position_last_geiger_hit[i].x()
						  << ' ' << collection_position_last_geiger_hit[i].y()
						  << ' ' << collection_position_last_geiger_hit[i].z() << std::endl;
				      }
				    calo_1.geiger_association = true;
				  }
			      }

			    if (calo_1.geiger_association)
			      {
				one_calo_one_track_energy_TH1F->Fill(calo_hit_map.begin()->second.energy * 1000);
				for (std::map<geomtools::geom_id, calo_hit_summary>::const_iterator it_calo = calo_hit_map.begin(); it_calo != calo_hit_map.end(); it_calo++)
				  {
				    int column = it_calo->first.get(2);
				    int row = it_calo->first.get(3);
				    one_calo_one_track_distribution_TH2F->Fill(column, row);
				  }

				for (std::set<geomtools::geom_id>::const_iterator it_geiger = geiger_hit_set.begin(); it_geiger != geiger_hit_set.end(); it_geiger++)
				  {
				    int layer = it_geiger->get(2);
				    int row   = it_geiger->get(3);
				    one_calo_one_track_geiger_distribution_TH2F->Fill(row, layer);
				  }
			      }
			    else
			      {
				if (is_display) std::clog << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& " << event_id << " &&&&&&&&&&&&&&&&" << std::endl;
			      }
			  }
		      } // end of number_of_main_calo == 1


		    // Two calos only :
		    if (number_of_main_calo == 2)
		      {
			calo_hit_summary calo_1 = calo_hit_map.begin()->second;
			calo_hit_summary calo_2 = calo_hit_map.rbegin()->second;
			if (is_display) std::clog << "Calo 1 GID : " << calo_1.geom_id << " Energy : " << calo_1.energy << " Time : " << calo_1.time << std::endl;
			if (is_display) std::clog << "Calo 2 GID : " << calo_2.geom_id << " Energy : " << calo_2.energy << " Time : " << calo_2.time << std::endl;
			double delta_time_2_calo = std::abs(calo_1.time - calo_2.time);
			double calo_energy_min = calo_1.energy;
			if (calo_energy_min > calo_2.energy) calo_energy_min = calo_2.energy;
			double calo_energy_max = calo_1.energy;
			if (calo_energy_max < calo_2.energy) calo_energy_max = calo_2.energy;
			double calo_energy_total = calo_energy_min + calo_energy_max;

			// 2 calos track or not
			two_calo_energy_min_TH1F->Fill(calo_energy_min * 1000);
			two_calo_energy_max_TH1F->Fill(calo_energy_max * 1000);
			two_calo_total_energy_TH1F->Fill(calo_energy_total * 1000);
			two_calo_delta_time_TH1F->Fill(delta_time_2_calo);

			for (std::map<geomtools::geom_id, calo_hit_summary>::const_iterator it_calo = calo_hit_map.begin(); it_calo != calo_hit_map.end(); it_calo++)
			  {
			    int column = it_calo->first.get(2);
			    int row = it_calo->first.get(3);
			    two_calo_distribution_TH2F->Fill(column, row);
			  }

			for (std::set<geomtools::geom_id>::const_iterator it_geiger = geiger_hit_set.begin(); it_geiger != geiger_hit_set.end(); it_geiger++)
			  {
			    int layer = it_geiger->get(2);
			    int row   = it_geiger->get(3);
			    two_calo_geiger_distribution_TH2F->Fill(row, layer);
			  }

			// 2 calos no track at all :
			if (number_of_layer == 0)
			  {
			    two_calo_no_track_energy_min_TH1F->Fill(calo_energy_min * 1000);
			    two_calo_no_track_energy_max_TH1F->Fill(calo_energy_max * 1000);
			    two_calo_no_track_total_energy_TH1F->Fill(calo_energy_total * 1000);
			    two_calo_no_track_delta_time_TH1F->Fill(delta_time_2_calo);

			    for (std::map<geomtools::geom_id, calo_hit_summary>::const_iterator it_calo = calo_hit_map.begin(); it_calo != calo_hit_map.end(); it_calo++)
			      {
				int column = it_calo->first.get(2);
				int row = it_calo->first.get(3);
				two_calo_no_track_distribution_TH2F->Fill(column, row);
			      }
			  }

			// 2 calos one track :
			if (full_track_event)
			  {
			    two_calos_one_full_track_writer.process(ER);
			    int side = 1;
			    int column_1 = calo_1.geom_id.get(2);
			    int row_1 = calo_1.geom_id.get(3);
			    int column_2 = calo_2.geom_id.get(2);
			    int row_2 = calo_2.geom_id.get(3);

			    double calo_1_zposition = calo_locator.get_row_z(side, row_1);
			    double calo_2_zposition = calo_locator.get_row_z(side, row_2);

			    for (unsigned int i = 0; i < collection_position_last_geiger_hit.size(); i++)
			      {
				if (calo_1_zposition > collection_position_last_geiger_hit[i].z() - 185
				    && calo_1_zposition < collection_position_last_geiger_hit[i].z() + 185)
				  {
				    if (is_display)
				      {
					std::clog << "Calorimeter and Geiger association" << std::endl;
					std::clog << "Calo 1 Z: " << calo_1_zposition << std::endl;
					std::clog << "GG hit #"<< i << ' ' << collection_position_last_geiger_hit[i].x()
						  << ' ' << collection_position_last_geiger_hit[i].y()
						  << ' ' << collection_position_last_geiger_hit[i].z() << std::endl;
				      }
				    calo_1.geiger_association = true;
				  }

				if (calo_2_zposition > collection_position_last_geiger_hit[i].z() - 185
				    && calo_2_zposition < collection_position_last_geiger_hit[i].z() + 185)
				  {
				    if (is_display)
				      {
					std::clog << "Calorimeter and Geiger association" << std::endl;
					std::clog << "Calo 2 Z: " << calo_2_zposition << std::endl;
					std::clog << "GG hit #"<< i << ' ' << collection_position_last_geiger_hit[i].x()
						  << ' ' << collection_position_last_geiger_hit[i].y()
						  << ' ' << collection_position_last_geiger_hit[i].z() << std::endl;
				      }
				    calo_2.geiger_association = true;
				  }
			      }

			    if (calo_1.geiger_association || calo_2.geiger_association)
			      {
				two_calo_one_track_energy_min_TH1F->Fill(calo_energy_min * 1000);
				two_calo_one_track_energy_max_TH1F->Fill(calo_energy_max * 1000);
				two_calo_one_track_total_energy_TH1F->Fill(calo_energy_total * 1000);
				two_calo_one_track_delta_time_TH1F->Fill(delta_time_2_calo);

				if (calo_1.geiger_association)
				  {
				    two_calo_one_track_electron_energy_TH1F->Fill(calo_1.energy * 1000);
				    two_calo_one_track_calo_interaction_distribution_TH2F->Fill(calo_1.left_most_hit_position.getX(), calo_1.left_most_hit_position.getZ());
				    if (is_display) std::clog << "Calo 1  X = " << calo_1.left_most_hit_position.getX()
							      << " Y = " << calo_1.left_most_hit_position.getY()
							      << " Z = " << calo_1.left_most_hit_position.getZ() << std::endl;

				  }
				else two_calo_one_track_gamma_energy_TH1F ->Fill(calo_1.energy * 1000);

				if (calo_2.geiger_association)
				  {
				    two_calo_one_track_electron_energy_TH1F->Fill(calo_2.energy * 1000);
				    two_calo_one_track_calo_interaction_distribution_TH2F->Fill(calo_2.left_most_hit_position.getX(), calo_2.left_most_hit_position.getZ());
				    if (is_display) std::clog << "Calo 2  X = " << calo_2.left_most_hit_position.getX()
							      << " Y = " << calo_2.left_most_hit_position.getY()
							      << " Z = " << calo_2.left_most_hit_position.getZ() << std::endl;
				  }
				else two_calo_one_track_gamma_energy_TH1F ->Fill(calo_2.energy * 1000);

				for (std::map<geomtools::geom_id, calo_hit_summary>::const_iterator it_calo = calo_hit_map.begin(); it_calo != calo_hit_map.end(); it_calo++)
				  {
				    int column = it_calo->first.get(2);
				    int row = it_calo->first.get(3);
				    two_calo_one_track_distribution_TH2F->Fill(column, row);
				  }

				for (std::set<geomtools::geom_id>::const_iterator it_geiger = geiger_hit_set.begin(); it_geiger != geiger_hit_set.end(); it_geiger++)
				  {
				    int layer = it_geiger->get(2);
				    int row   = it_geiger->get(3);
				    two_calo_one_track_geiger_distribution_TH2F->Fill(row, layer);
				  }

				// Calcul the angular between source vertex and middle calorimeter hits :
				geomtools::vector_3d vertex_position = flaged_SD.get_vertex();
				geomtools::vector_3d calo_1_position = calo_locator.get_block_position(side, column_1, row_1);
				geomtools::vector_3d calo_2_position = calo_locator.get_block_position(side, column_2, row_2);

				double angle_calo_1 = (180.0 / 3.14159265358979323846) * vertex_position.angle(calo_1_position);
				double angle_calo_2 = (180.0 / 3.14159265358979323846) * vertex_position.angle(calo_2_position);
				two_calo_one_track_angle_TH1F->Fill(angle_calo_1);
				two_calo_one_track_angle_TH1F->Fill(angle_calo_2);

			      } // end of calo geiger association
			  } // end of is full track
		      } // end of 2 calos

		  } // end of number of calo > 0
	      } // end of match rules

	  } // end of if ER has flaged_SD_bank_label

	event_id++;

	ER.clear();
      } // end of reader is terminated

    root_file->cd();

    root_file->Write();
    root_file->Close();

    std::clog << "The end." << std::endl;
  } // end of try

  catch (std::exception & error) {
    DT_LOG_FATAL(logging, error.what());
    error_code = EXIT_FAILURE;
  }

  catch (...) {
    DT_LOG_FATAL(logging, "Unexpected error!");
    error_code = EXIT_FAILURE;
  }

  falaise::terminate();
  return error_code;
}
