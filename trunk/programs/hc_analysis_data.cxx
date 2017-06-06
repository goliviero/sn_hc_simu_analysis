// hc_analysis_data.cxx
// Standard libraries :
#include <bitset>
// #include <iostream>

// Third party:
// - Boost:
#include <boost/program_options.hpp>

// - Bayeux/datatools:
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

// Root :
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

// This project :
#include "data_statistics_simu.hpp"
#include "hc_constants.hpp"

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

  try {

    std::vector<std::string> input_filenames; // = "";
    std::string output_path = "";
    std::string calo_mapping_config = "";
    std::string tracker_mapping_config = "";
    std::size_t max_events  = 0;
    bool        is_debug    = false;
    double      calo_threshold_kev  = 0;

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
      ("calo-threshold,c",
       po::value<double>(& calo_threshold_kev)->default_value(hc_constants::CALO_COMMISSIONING_HIGH_THRESHOLD_KEV),
       "set the calorimeter threshold in keV")
      ("calo_mapping,C",
       po::value<std::string>(& calo_mapping_config),
       "set the calorimeter mapping configuration from a datatools::properties ASCII file")
      ("tracker_mapping,T",
       po::value<std::string>(& tracker_mapping_config),
       "set the tracker mapping configuration from a datatools::properties ASCII file")
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
      return(1);
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

    std::clog << "INFO : Program sorting raw data file to reduce the number of event / file and produce root histograms for analysis simulated data from half commissioning" << std::endl;

    std::string manager_config_file;
    manager_config_file = "@falaise:config/snemo/demonstrator/geometry/4.0/manager.conf";
    datatools::fetch_path_with_env(manager_config_file);
    datatools::properties manager_config;
    datatools::properties::read_config (manager_config_file,
					manager_config);
    geomtools::manager my_geom_manager;
    manager_config.update ("build_mapping", true);
    if (manager_config.has_key ("mapping.excluded_categories"))	manager_config.erase ("mapping.excluded_categories");
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

    // Calo and tracker half commissioning mapping configuration :
    datatools::properties calo_config;
    if (!calo_mapping_config.empty()) {
      calo_config.read_configuration(calo_mapping_config);
    }

    datatools::properties tracker_config;
    if (!tracker_mapping_config.empty()) {
      tracker_config.read_configuration(tracker_mapping_config);
    }

    geomtools::id_selector hc_calo_selector(my_geom_manager.get_id_mgr());
    if (!calo_mapping_config.empty()) {
      std::ifstream ifile(tracker_mapping_config);
      bool empty = false;
      if (ifile.peek() == std::ifstream::traits_type::eof()) empty = true;
      if (!empty)  hc_calo_selector.initialize(calo_config);
    }
    if (is_debug) hc_calo_selector.dump(std::clog, "Half commissioning calo selector: ");

    geomtools::id_selector hc_geiger_selector(my_geom_manager.get_id_mgr());
    if (!tracker_mapping_config.empty()) {
      std::ifstream ifile(tracker_mapping_config);
      bool empty = false;
      if (ifile.peek() == std::ifstream::traits_type::eof()) empty = true;
      if (!empty) hc_geiger_selector.initialize(tracker_config);
    }
    if (is_debug) hc_geiger_selector.dump(std::clog, "Half commissioning Geiger selector: ");

    //==============================================//
    //          output files  and writers           //
    //==============================================//

    // Name of full track (9 layers hit) SD output file :
    std::string calo_tracker_events_brio = output_path + "output_calo_tracker_events.brio";

    // Event writer for full track (9 layers hit) :
    dpp::output_module calo_tracker_events_writer;
    datatools::properties calo_tracker_events_config;
    calo_tracker_events_config.store ("logging.priority", "fatal");
    calo_tracker_events_config.store ("files.mode", "single");
    calo_tracker_events_config.store ("files.single.filename", calo_tracker_events_brio);
    calo_tracker_events_writer.grab_metadata_store() = iMetadataStore;
    calo_tracker_events_writer.initialize_standalone(calo_tracker_events_config);

    // Output ROOT file :
    std::string string_buffer = output_path + "output_rootfile.root";
    datatools::fetch_path_with_env(string_buffer);

    TFile* root_file = new TFile(string_buffer.c_str(), "RECREATE");

    // Event counter :
    int event_id    = 0;

    data_statistics_simu my_dss;
    my_dss.initialize();

    while (!reader.is_terminated())
      {
	DT_LOG_DEBUG(logging, "Event #" << event_id);
	reader.process(ER);

	bool full_track_event = false;

	// A plain `mctools::simulated_data' object is stored here :
	if (ER.has(SD_bank_label) && ER.is_a<mctools::simulated_data>(SD_bank_label))
	  {
	    // Access to the "SD" bank with a stored `mctools::simulated_data' :
	    const mctools::simulated_data & SD = ER.get<mctools::simulated_data>(SD_bank_label);

	    // New SD bank flaged for Geiger cell already hit
	    mctools::simulated_data flaged_SD = SD;

	    // First loop on all hits to merge each calo hit in the same OM, also tag GG cells hit several times :

	    std::map<geomtools::geom_id, calo_hit_summary> calo_hit_map;
	    // If main calo hits :
	    if (flaged_SD.has_step_hits("calo"))
	      {
		mctools::simulated_data::hit_handle_collection_type BSHC = flaged_SD.get_step_hits("calo");
		DT_LOG_DEBUG(logging, "BSCH calo step hits # = " << BSHC.size());

		for (mctools::simulated_data::hit_handle_collection_type::const_iterator i = BSHC.begin();
		     i != BSHC.end();
		     i++)
		  {
		    const mctools::base_step_hit & BSH = i->get();
		    // extract the corresponding geom ID:
		    const geomtools::geom_id & main_calo_gid = BSH.get_geom_id();

		    if (hc_calo_selector.match(main_calo_gid)) {

		      bool calo_hit_is_in_map = calo_hit_map.find(main_calo_gid) != calo_hit_map.end();
		      //int column = main_calo_gid.get(2);
		      // int row = main_calo_gid.get(3);

		      if (!calo_hit_is_in_map)
			{
			  // Create a new calo hit and insert in the map only if it match commissioning rules
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

		    } // end of match rules selection

		  } // end of for i BSHC

		// Check if all calo summary hit pass the threshold :
		for (std::map<geomtools::geom_id, calo_hit_summary>::iterator it_calo = calo_hit_map.begin();
		     it_calo != calo_hit_map.end();)
		  {
		    if (it_calo->second.energy * 1000 < calo_threshold_kev) it_calo = calo_hit_map.erase(it_calo);
		    else it_calo++;
		  }

	      } // end of if has step hits "calo"

	    std::set<geomtools::geom_id> geiger_hit_set;
	    std::vector<geomtools::vector_3d> collection_position_last_geiger_hit;
	    std::size_t geiger_last_layer = hc_constants::NUMBER_OF_GEIGER_LAYERS - 1;

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
		DT_LOG_DEBUG(logging, "BSCH geiger step hits # = " << BSHC_gg.size());
		for (mctools::simulated_data::hit_handle_collection_type::const_iterator i = BSHC_gg.begin();
		     i != BSHC_gg.end();
		     i++)
		  {
		    const mctools::base_step_hit & BSH = i->get();
		    // if (is_debug) BSH.tree_dump(std::clog, "A Geiger Base Step Hit : ", "INFO : ");

		    // Second loop on Geiger hits, ignore cells hit 2 or more times.
		    if (BSH.get_auxiliaries().has_flag("geiger_already_hit") || BSH.get_auxiliaries().has_flag("other_geiger_already_hit")) {}
		    else
		      {
			const geomtools::geom_id & geiger_gid = BSH.get_geom_id();

			geiger_hit_set.insert(geiger_gid);
			// If the last Geiger at layer 8 is hit, push back position for calorimeter association
			if (geiger_gid.get(2) == geiger_last_layer) // last layer
			  {
			    collection_position_last_geiger_hit.push_back(BSH.get_position_stop());
			  }
		      }
		  } // end of for
	      } // end of if has step hits "gg"

	    /*******************************************************/
	    /* Begining analysis (based on calo map and geiger set */
	    /*******************************************************/

	    // Fill histograms in ROOT file :

	    // int number_of_main_calo = calo_hit_map.size();
	    // int number_of_geiger    = geiger_hit_set.size();

	    // For each calorimeter, add it in the histogram
	    // Calorimeter 'exists' only if E_calo > threshold
	    double calo_tref;
	    datatools::invalidate(calo_tref);

	    double total_energy = 0;
	    for (std::map<geomtools::geom_id, calo_hit_summary>::const_iterator it_calo = calo_hit_map.begin();
		 it_calo != calo_hit_map.end();
		 it_calo++)
	      {
		if (it_calo == calo_hit_map.begin()) calo_tref = it_calo->second.time;
		if (it_calo->second.time < calo_tref) calo_tref = it_calo->second.time;

		int column = it_calo->first.get(2);
		int row = it_calo->first.get(3);
		my_dss.calo_distrib_ht_TH2F->Fill(column, row);
		my_dss.calo_ht_energy_TH1F[column]->Fill(it_calo->second.energy * 1000);
		total_energy+=it_calo->second.energy;
	      }

	    my_dss.calo_ht_total_energy_TH1F->Fill(total_energy * 1000);

	    std::bitset<hc_constants::NUMBER_OF_GEIGER_LAYERS> layer_projection = 0x0;

	    // For each Geiger cell, add it in the histogram
	    for (std::set<geomtools::geom_id>::const_iterator it_geiger = geiger_hit_set.begin();
		 it_geiger != geiger_hit_set.end();
		 it_geiger++)
	      {
		int layer = it_geiger->get(2);
		int row   = it_geiger->get(3);
		// if (is_debug) std::clog << "Layer  = " << layer << " Row = " << row << std::endl;
		my_dss.tracker_total_distribution_TH2F->Fill(row, layer);
		layer_projection.set(layer, true);
	      }
	    int number_of_layer = layer_projection.count();
	    if (number_of_layer == hc_constants::NUMBER_OF_GEIGER_LAYERS) full_track_event = true;

	    if (full_track_event) DT_LOG_DEBUG(logging, "Full track event !");

	    // // Two calos only :
	    // if (number_of_main_calo == 2)
	    //   {
	    // 	calo_hit_summary calo_1 = calo_hit_map.begin()->second;
	    // 	calo_hit_summary calo_2 = calo_hit_map.rbegin()->second;
	    // 	if (is_debug) std::clog << "Calo 1 GID : " << calo_1.geom_id << " Energy : " << calo_1.energy << " Time : " << calo_1.time << std::endl;
	    // 	if (is_debug) std::clog << "Calo 2 GID : " << calo_2.geom_id << " Energy : " << calo_2.energy << " Time : " << calo_2.time << std::endl;
	    // 	double delta_time_2_calo = std::abs(calo_1.time - calo_2.time);
	    // 	double calo_energy_min = calo_1.energy;
	    // 	if (calo_energy_min > calo_2.energy) calo_energy_min = calo_2.energy;
	    // 	double calo_energy_max = calo_1.energy;
	    // 	if (calo_energy_max < calo_2.energy) calo_energy_max = calo_2.energy;
	    // 	double calo_energy_total = calo_energy_min + calo_energy_max;

	    // 	for (std::map<geomtools::geom_id, calo_hit_summary>::const_iterator it_calo = calo_hit_map.begin(); it_calo != calo_hit_map.end(); it_calo++)
	    // 	  {
	    // 	    int column = it_calo->first.get(2);
	    // 	    int row = it_calo->first.get(3);
	    // 	    two_calo_distribution_TH2F->Fill(column, row);
	    // 	  }

	    // 	for (std::set<geomtools::geom_id>::const_iterator it_geiger = geiger_hit_set.begin(); it_geiger != geiger_hit_set.end(); it_geiger++)
	    // 	  {
	    // 	    int layer = it_geiger->get(2);
	    // 	    int row   = it_geiger->get(3);
	    // 	    two_calo_geiger_distribution_TH2F->Fill(row, layer);
	    // 	  }

	    // 	// 2 calos no track at all :
	    // 	if (number_of_layer == 0)
	    // 	  {
	    // 	    two_calo_no_track_energy_min_TH1F->Fill(calo_energy_min * 1000);
	    // 	    two_calo_no_track_energy_max_TH1F->Fill(calo_energy_max * 1000);
	    // 	    two_calo_no_track_total_energy_TH1F->Fill(calo_energy_total * 1000);
	    // 	    two_calo_no_track_delta_time_TH1F->Fill(delta_time_2_calo);

	    // 	    for (std::map<geomtools::geom_id, calo_hit_summary>::const_iterator it_calo = calo_hit_map.begin(); it_calo != calo_hit_map.end(); it_calo++)
	    // 	      {
	    // 		int column = it_calo->first.get(2);
	    // 		int row = it_calo->first.get(3);
	    // 		two_calo_no_track_distribution_TH2F->Fill(column, row);
	    // 	      }
	    // 	  }

	    // 	// 2 calos one track :
	    // 	if (full_track_event)
	    // 	  {
	    // 	    int side = 1;
	    // 	    int column_1 = calo_1.geom_id.get(2);
	    // 	    int row_1 = calo_1.geom_id.get(3);
	    // 	    int column_2 = calo_2.geom_id.get(2);
	    // 	    int row_2 = calo_2.geom_id.get(3);

	    // 	    double calo_1_zposition = calo_locator.get_row_z(side, row_1);
	    // 	    double calo_2_zposition = calo_locator.get_row_z(side, row_2);

	    // 	    for (unsigned int i = 0; i < collection_position_last_geiger_hit.size(); i++)
	    // 	      {
	    // 		if (calo_1_zposition > collection_position_last_geiger_hit[i].z() - 185
	    // 		    && calo_1_zposition < collection_position_last_geiger_hit[i].z() + 185)
	    // 		  {
	    // 		    if (is_debug)
	    // 		      {
	    // 			std::clog << "Calorimeter and Geiger association" << std::endl;
	    // 			std::clog << "Calo 1 Z: " << calo_1_zposition << std::endl;
	    // 			std::clog << "GG hit #"<< i << ' ' << collection_position_last_geiger_hit[i].x()
	    // 				  << ' ' << collection_position_last_geiger_hit[i].y()
	    // 				  << ' ' << collection_position_last_geiger_hit[i].z() << std::endl;
	    // 		      }
	    // 		    calo_1.geiger_association = true;
	    // 		  }

	    // 		if (calo_2_zposition > collection_position_last_geiger_hit[i].z() - 185
	    // 		    && calo_2_zposition < collection_position_last_geiger_hit[i].z() + 185)
	    // 		  {
	    // 		    if (is_debug)
	    // 		      {
	    // 			std::clog << "Calorimeter and Geiger association" << std::endl;
	    // 			std::clog << "Calo 2 Z: " << calo_2_zposition << std::endl;
	    // 			std::clog << "GG hit #"<< i << ' ' << collection_position_last_geiger_hit[i].x()
	    // 				  << ' ' << collection_position_last_geiger_hit[i].y()
	    // 				  << ' ' << collection_position_last_geiger_hit[i].z() << std::endl;
	    // 		      }
	    // 		    calo_2.geiger_association = true;
	    // 		  }
	    // 	      }

	    // 	    if (calo_1.geiger_association || calo_2.geiger_association)
	    // 	      {
	    // 		two_calo_one_track_energy_min_TH1F->Fill(calo_energy_min * 1000);
	    // 		two_calo_one_track_energy_max_TH1F->Fill(calo_energy_max * 1000);
	    // 		two_calo_one_track_total_energy_TH1F->Fill(calo_energy_total * 1000);
	    // 		two_calo_one_track_delta_time_TH1F->Fill(delta_time_2_calo);

	    // 		if (calo_1.geiger_association)
	    // 		  {
	    // 		    two_calo_one_track_electron_energy_TH1F->Fill(calo_1.energy * 1000);
	    // 		    two_calo_one_track_calo_interaction_distribution_TH2F->Fill(calo_1.left_most_hit_position.getX(), calo_1.left_most_hit_position.getZ());
	    // 		    if (is_debug) std::clog << "Calo 1  X = " << calo_1.left_most_hit_position.getX()
	    // 					    << " Y = " << calo_1.left_most_hit_position.getY()
	    // 					    << " Z = " << calo_1.left_most_hit_position.getZ() << std::endl;

	    // 		  }
	    // 		else two_calo_one_track_gamma_energy_TH1F ->Fill(calo_1.energy * 1000);

	    // 		if (calo_2.geiger_association)
	    // 		  {
	    // 		    two_calo_one_track_electron_energy_TH1F->Fill(calo_2.energy * 1000);
	    // 		    two_calo_one_track_calo_interaction_distribution_TH2F->Fill(calo_2.left_most_hit_position.getX(), calo_2.left_most_hit_position.getZ());
	    // 		    if (is_debug) std::clog << "Calo 2  X = " << calo_2.left_most_hit_position.getX()
	    // 					    << " Y = " << calo_2.left_most_hit_position.getY()
	    // 					    << " Z = " << calo_2.left_most_hit_position.getZ() << std::endl;
	    // 		  }
	    // 		else two_calo_one_track_gamma_energy_TH1F ->Fill(calo_2.energy * 1000);

	    // 		for (std::map<geomtools::geom_id, calo_hit_summary>::const_iterator it_calo = calo_hit_map.begin(); it_calo != calo_hit_map.end(); it_calo++)
	    // 		  {
	    // 		    int column = it_calo->first.get(2);
	    // 		    int row = it_calo->first.get(3);
	    // 		    two_calo_one_track_distribution_TH2F->Fill(column, row);
	    // 		  }

	    // 		for (std::set<geomtools::geom_id>::const_iterator it_geiger = geiger_hit_set.begin(); it_geiger != geiger_hit_set.end(); it_geiger++)
	    // 		  {
	    // 		    int layer = it_geiger->get(2);
	    // 		    int row   = it_geiger->get(3);
	    // 		    two_calo_one_track_geiger_distribution_TH2F->Fill(row, layer);
	    // 		  }

	    // 		// Calcul the angular between source vertex and middle calorimeter hits :
	    // 		geomtools::vector_3d vertex_position = flaged_SD.get_vertex();
	    // 		geomtools::vector_3d calo_1_position = calo_locator.get_block_position(side, column_1, row_1);
	    // 		geomtools::vector_3d calo_2_position = calo_locator.get_block_position(side, column_2, row_2);

	    // 		double angle_calo_1 = (180.0 / 3.14159265358979323846) * vertex_position.angle(calo_1_position);
	    // 		double angle_calo_2 = (180.0 / 3.14159265358979323846) * vertex_position.angle(calo_2_position);
	    // 		two_calo_one_track_angle_TH1F->Fill(angle_calo_1);
	    // 		two_calo_one_track_angle_TH1F->Fill(angle_calo_2);

	    // 	      } // end of calo geiger association
	    // 	  } // end of is full track
	    //   } // end of 2 calos

	  } // end of if ER has flaged_SD_bank_label

	event_id++;

	ER.clear();
      } // end of reader is terminated

    root_file->cd();
    my_dss.save_in_root_file(root_file);
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
