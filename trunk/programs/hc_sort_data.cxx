// hc_sort_data.cxx
// Standard libraries :
// #include <iostream>
// #include <bitset>
// #include <fstream>

// Third party:
// - Boost:
#include <boost/program_options.hpp>

// // - Bayeux/datatools:
#include <datatools/utils.h>
#include <datatools/io_factory.h>

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

int main( int  argc_ , char **argv_  )
{
  falaise::initialize(argc_, argv_);
  int error_code = EXIT_SUCCESS;
  datatools::logger::priority logging = datatools::logger::PRIO_FATAL;

  try {

    std::vector<std::string> input_filenames; // = "";
    std::string output_path = "";
    std::size_t max_events  = 0;
    bool is_debug = false;

    // Parse options:
    namespace po = boost::program_options;
    po::options_description opts("Allowed options");
    opts.add_options()
      ("help,h", "produce help message")
      ("debug,d", "debug mode")
      ("input,i",
       po::value<std::vector<std::string> >(& input_filenames)->multitoken(),
       "set a list of input files")
      ("output,o",
       po::value<std::string>(& output_path),
       "set the output path")
      ("number_events,n",
       po::value<std::size_t>(& max_events)->default_value(10),
       "set the maximum number of events")
      ; // end of options description

    // Describe command line arguments :
    po::variables_map vm;
    po::store(po::command_line_parser(argc_, argv_)
	      .options(opts)
	      .run(), vm);
    po::notify(vm);

    // Use command line arguments :
    if (vm.count("help")) {
      std::cout << "Usage : " << std::endl;
      std::cout << opts << std::endl;

      return(error_code);
    }

    // Use command line arguments :
    else if (vm.count("debug")) {
     is_debug = true;
    }

    if (is_debug) logging = datatools::logger::PRIO_DEBUG;

    DT_LOG_INFORMATION(logging, "List of input file(s) : ");
    for (auto file = input_filenames.begin();
	 file != input_filenames.end();
	 file++) std::clog << *file << ' ';
    DT_THROW_IF(input_filenames.size() == 0, std::logic_error, "No input file(s) ! ");

    DT_LOG_INFORMATION(logging, "Output path for files = " + output_path);
    if (output_path.empty()) {
      output_path = ".";
      DT_LOG_INFORMATION(logging, "No output path, default output path is = " + output_path);
    }

    std::clog << "INFO : Welcome in the Half Commissioning sorting data program" << std::endl;
    std::clog << "INFO : Program sorting raw data file to reduce the number of event / file" << std::endl;
    std::clog << "INFO : Events have to touch GG cells or OMs in given zones to be saved" << std::endl;

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

    // Simulated Data "SD" bank label :
    std::string SD_bank_label = "SD";

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

    int max_record_total = static_cast<int>(max_events) * static_cast<int>(input_filenames.size());
    std::clog << "max_record total = " << max_record_total << std::endl;
    std::clog << "max_events       = " << max_events << std::endl;


    // Event reader :
    dpp::input_module reader;
    datatools::properties reader_config;
    reader_config.store ("logging.priority", "debug");
    reader_config.store("files.mode", "list");
    reader_config.store("files.list.filenames", input_filenames);
    reader_config.store("max_record_total", max_record_total);
    reader_config.store("max_record_per_file", static_cast<int>(max_events));
    reader.initialize_standalone (reader_config);
    datatools::multi_properties iMetadataStore = reader.get_metadata_store();
    // reader.tree_dump(std::clog, "Simulated data reader module");

    // Event record :
    datatools::things ER;

    /********************************************************************/
    /* DESIGN OF THE SOFT : JUST FOR ZONE 2 FOR THE MOMENT (no need to  */
    /* simulate all half zones (unefficient) maybe later. TO DO : adapt */
    /* the software for all half zones (take into account row/column or */
    /* free spots vertexes (row / column approach needs sometimes to be */
    /* process 2 times (column 3 4 5 count for zone 0 and 1 for example */
    /********************************************************************/

    // unsigned int hc_half_zone_number = 0;

    // Selection relative to the commissioning mapping (calo + tracker)
    std::size_t hc_calo_column = 1; // column calo was 1 (real condition : PMT 0 and 12 inactive)
    std::size_t hc_geiger_row = 4;  // row was 4 (real condition : cell 1 and 8 inactive)

    // Id selector rules for half commissioning mapping :
    std::string hc_calo_rules = "category='calorimeter_block' module={0} side={1} column={" + std::to_string(hc_calo_column) + "} row={*} part={*}";
    // Only column 4 is sensitive because commissioning
    std::string hc_geiger_rules = "category='drift_cell_core' module={0} side={1} layer={*} row={"
      + std::to_string(hc_geiger_row) + "}";

    geomtools::id_selector hc_calo_selector(my_geom_manager.get_id_mgr());
    hc_calo_selector.initialize(hc_calo_rules);
    if (is_debug) hc_calo_selector.dump(std::clog, "Half commissioning calo selector: ");

    geomtools::id_selector hc_geiger_selector(my_geom_manager.get_id_mgr());
    hc_geiger_selector.initialize(hc_geiger_rules);
    if (is_debug) hc_geiger_selector.dump(std::clog, "Half commissioning Geiger selector: ");

    //============================================//
    //          output file  and writer           //
    //============================================//

    // Name of sorted (matching rules) SD output file :
    std::string sorted_sd_brio = output_path + "match_rules.brio";

    // Event writer for sorted SD :
    dpp::output_module sorted_writer;
    datatools::properties sorted_writer_config;
    sorted_writer_config.store ("logging.priority", "fatal");
    sorted_writer_config.store ("files.mode", "single");
    sorted_writer_config.store ("files.single.filename", sorted_sd_brio);
    sorted_writer.grab_metadata_store() = iMetadataStore;
    sorted_writer.initialize_standalone(sorted_writer_config);

    // Name of sorted (matching rules) SD output file :
    std::string sorted_with_geiger_brio = output_path + "match_rules_with_geiger.brio";

    // Event writer for sorted SD :
    dpp::output_module sorted_with_geiger_writer;
    datatools::properties sorted_with_geiger_config;
    sorted_with_geiger_config.store ("logging.priority", "fatal");
    sorted_with_geiger_config.store ("files.mode", "single");
    sorted_with_geiger_config.store ("files.single.filename", sorted_with_geiger_brio);
    sorted_with_geiger_writer.grab_metadata_store() = iMetadataStore;
    sorted_with_geiger_writer.initialize_standalone(sorted_with_geiger_config);

    // Event counter :
    int event_id    = 0;

    while (!reader.is_terminated())
      {
	DT_LOG_DEBUG(logging, "Event #" + event_id);
	reader.process(ER);

	bool match_rules_event = false;
	bool match_rules_with_geiger = false;

	// A plain `mctools::simulated_data' object is stored here :
	if (ER.has(SD_bank_label) && ER.is_a<mctools::simulated_data>(SD_bank_label))
	  {
	    // Access to the "SD" bank with a stored `mctools::simulated_data' :
	    const mctools::simulated_data & SD = ER.get<mctools::simulated_data>(SD_bank_label);

	    // Main calo hits :
	    if (SD.has_step_hits("calo"))
	      {
		mctools::simulated_data::hit_handle_collection_type BSHC = SD.get_step_hits("calo");
		if (is_debug) std::clog << "BSCH calo step hits # = " << BSHC.size() << std::endl;

		for (mctools::simulated_data::hit_handle_collection_type::const_iterator i = BSHC.begin();
		     i != BSHC.end();
		     i++)
		  {
		    const mctools::base_step_hit & BSH = i->get();
		    // extract the corresponding geom ID:
		    const geomtools::geom_id & main_calo_gid = BSH.get_geom_id();
		    if (hc_calo_selector.match(main_calo_gid)) match_rules_event = true;

		  } // end of for i BSHC

	      } // end of if has step hits "calo"


	    if (SD.has_step_hits("gg"))
	      {
		mctools::simulated_data::hit_handle_collection_type BSHC_gg = SD.get_step_hits("gg");
		if (is_debug) std::clog << "BSCH geiger step hits # = " << BSHC_gg.size() << std::endl;
		for (mctools::simulated_data::hit_handle_collection_type::const_iterator i = BSHC_gg.begin();
		     i != BSHC_gg.end();
		     i++)
		  {
		    const mctools::base_step_hit & BSH = i->get();
		    if (is_debug) BSH.tree_dump(std::clog, "A Geiger Base Step Hit : ", "INFO : ");
		    const geomtools::geom_id & geiger_gid = BSH.get_geom_id();

		    if (hc_geiger_selector.match(geiger_gid))
		      {
			match_rules_event = true;
			match_rules_with_geiger = true;
		      }

		  } // end of i bsh
	      } // end of if has step hits "gg"

	    if (match_rules_event) sorted_writer.process(ER);
	    if (match_rules_with_geiger) sorted_with_geiger_writer.process(ER);

	  } // end of if ER has flaged_SD_bank_label

	event_id++;

	ER.clear();
      } // end of reader is terminated


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
