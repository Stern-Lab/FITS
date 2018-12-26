/*
    FITS - Flexible Inference from Time-Series data
    (c) 2016-2018 by Tal Zinger
    tal.zinger@outlook.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "ResultsStats.hpp"


void ResultsStats::CalculateStatsFitness()
{
    auto wt_allele = _result_vector[0].wt_index;
    _num_alleles = _result_vector[0].fitness_values.size();
    
    _num_results = _result_vector.size();
    
    _allele_Nu.resize(_num_alleles, 0.0f);
    levenes_pval.resize(_num_alleles, -1.0f);
    
    MATRIX_TYPE posterior_matrix( _num_results, _num_alleles );
    MATRIX_TYPE prior_matrix( _num_results, _num_alleles );
    
    
    
    boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > acc_distance;
    
    std::vector< boost::accumulators::accumulator_set<
    FLOAT_TYPE,
    boost::accumulators::stats<
    boost::accumulators::tag::variance,
    boost::accumulators::tag::mean,
    boost::accumulators::tag::median,
    boost::accumulators::tag::min,
    boost::accumulators::tag::max> > > acc_vec_fitness(_num_alleles);

    allele_MAD.resize(_num_alleles, 0.0f);
    allele_max_fitness.resize(_num_alleles, 0.0f);
    allele_mean_fitness.resize(_num_alleles, 0.0f);
    allele_sd_fitness.resize(_num_alleles, 0.0f);
    allele_median_fitness.resize(_num_alleles, 0.0f);
    allele_min_fitness.resize(_num_alleles, 0.0f);
    
    allele_min_95percentile_fitness.resize(_num_alleles, 0.0f);
    allele_max_95percentile_fitness.resize(_num_alleles, 0.0f);
    
    allele_pval.resize(_num_alleles, 0.0f);
    lethal_counter.resize(_num_alleles, 0.0f);
    deleterious_counter.resize(_num_alleles, 0.0f);
    neutral_counter.resize(_num_alleles, 0.0f);
    advantageous_counter.resize(_num_alleles, 0.0f);
    //unassigned_counter(_num_alleles, 0.0f);
    
    lethal_percent.resize(_num_alleles, 0.0f);
    deleterious_percent.resize(_num_alleles, 0.0f);
    neutral_percent.resize(_num_alleles, 0.0f);
    advantageous_percent.resize(_num_alleles, 0.0f);
    unassigned_percent.resize(_num_alleles, 0.0f);
    
    allele_category.resize(_num_alleles, AlleleCategory::Undefined);
    
    
    // measure of relaiability. we prefer Nu>=1
    // u = f(wt->current_allele)
    for (auto current_allele = 0; current_allele<_num_alleles; ++current_allele) {
        
        if ( _zparams.GetInt( "Debug", 0 ) > 0 ) {
            std::cout << "Calculating Nu. N=" << _result_vector[0].N <<
            "; u=" << _result_vector[0].mutation_rates(_result_vector[0].wt_index, current_allele) << std::endl;
        }
        
        _allele_Nu[current_allele] =
            _result_vector[0].N *
                _result_vector[0].mutation_rates(_result_vector[0].wt_index, current_allele);
    }
    
    
    
    std::vector< std::vector<FLOAT_TYPE> > allele_fitness_storage( _num_alleles );
    
    for (auto current_sim_index=0; current_sim_index<_result_vector.size(); ++current_sim_index ) {
        
        auto current_sim_data = _result_vector[current_sim_index];
        
        acc_distance(current_sim_data.distance_from_actual);
        
        for (auto current_allele = 0; current_allele<current_sim_data.fitness_values.size(); current_allele++) {
         
            allele_fitness_storage[current_allele].push_back( _result_vector[current_sim_index].fitness_values[current_allele] );
            
            auto tmpval = current_sim_data.fitness_values[current_allele];
            
            posterior_matrix( current_sim_index, current_allele ) = tmpval;
            
            acc_vec_fitness[current_allele](tmpval);
            
            
           
            // 2016-10-04 - removed categories for lethal and neutral
            if ( tmpval > FITNESS_NEUTRAL ) {
                ++advantageous_counter[current_allele];
            }
            else if ( tmpval < FITNESS_NEUTRAL ) {
                ++deleterious_counter[current_allele];
            }
            else {
                ++neutral_counter[current_allele];
            }
        } // for with _counters
        
    }
    
    
    
    
    _distance_min = boost::accumulators::min(acc_distance);
    _distance_max = boost::accumulators::max(acc_distance);
    _distance_mean = boost::accumulators::mean(acc_distance);
    _distance_sd = std::sqrt(boost::accumulators::variance(acc_distance));
    
    for (auto current_allele = 0; current_allele<_num_alleles; ++current_allele) {
        
        // added boost::accumulators to prevent any ambiguity with std library
        allele_mean_fitness[current_allele] = boost::accumulators::mean(acc_vec_fitness[current_allele]);
        allele_sd_fitness[current_allele] = sqrt(boost::accumulators::variance(acc_vec_fitness[current_allele]));
        allele_max_fitness[current_allele] = boost::accumulators::max(acc_vec_fitness[current_allele]);
        allele_min_fitness[current_allele] = boost::accumulators::min(acc_vec_fitness[current_allele]);
        
        
        // allele_median_fitness[current_allele] = GetMedian( allele_fitness_storage[current_allele] );
        allele_median_fitness[current_allele] = boost::accumulators::median(acc_vec_fitness[current_allele]);
        
        // calculate Bayesian pval - P(w<1|data)
        allele_pval[current_allele] =
        static_cast<double>(lethal_counter[current_allele] +
                            deleterious_counter[current_allele]) * 100.0 /
        static_cast<double>(_num_results);
    }
    
    
    // Calculate MAD for the alleles
    
    // start with distances
    //std::vector< std::vector<FLOAT_TYPE> > allele_distance_vec( _num_alleles );
    
    for ( auto current_allele=0; current_allele<allele_median_fitness.size(); ++current_allele ) {
        
        boost::accumulators::accumulator_set<
        FLOAT_TYPE,
        boost::accumulators::stats<
        boost::accumulators::tag::median > > acc_distance_for_mad;
        
        for ( auto current_fitness : allele_fitness_storage[current_allele] ) {
            // allele_distance_vec[current_allele].push_back( std::fabs(current_fitness - allele_median_fitness[current_allele] ) );
            // allele_MAD[current_allele] = GetMedian( allele_distance_vec[current_allele] );
            acc_distance_for_mad( std::fabs(current_fitness - allele_median_fitness[current_allele] ) );
        }
        
        allele_MAD[current_allele] = boost::accumulators::median(acc_distance_for_mad);
    }

    
    for (auto current_allele = 0; current_allele < _num_alleles; ++current_allele) {
        
        // first so it will be the default
        allele_category[current_allele] = AlleleCategory::Undefined;
        /*
         lethal_percent[current_allele] = lethal_counter[current_allele] * 100 / result_vector.size();
         if (lethal_percent[current_allele] >= CATEGORY_INCLUSION_STRICT_THRESHOLD)
         allele_category[current_allele] = AlleleCategory::Lethal;
         else if (lethal_percent[current_allele] >= CATEGORY_INCLUSION_RELAXED_THRESHOLD)
         allele_category[current_allele] = AlleleCategory::Possible_lethal;
         */
        deleterious_percent[current_allele] = deleterious_counter[current_allele] * 100.0f / static_cast<double>(_result_vector.size());
        if (deleterious_percent[current_allele] >= CATEGORY_INCLUSION_STRICT_THRESHOLD)
            allele_category[current_allele] = AlleleCategory::Deleterious;
        else if (deleterious_percent[current_allele] >= CATEGORY_INCLUSION_RELAXED_THRESHOLD)
            allele_category[current_allele] = AlleleCategory::Possible_deleterious;
        
        neutral_percent[current_allele] = neutral_counter[current_allele] * 100.0f / static_cast<double>(_result_vector.size());
        if (neutral_percent[current_allele] >= CATEGORY_INCLUSION_STRICT_THRESHOLD)
            allele_category[current_allele] = AlleleCategory::Neutral;
        else if (neutral_percent[current_allele] >= CATEGORY_INCLUSION_RELAXED_THRESHOLD)
            allele_category[current_allele] = AlleleCategory::Possible_neutral;
        
        advantageous_percent[current_allele] = advantageous_counter[current_allele] * 100.0f / static_cast<double>(_result_vector.size());
        if (advantageous_percent[current_allele] >= CATEGORY_INCLUSION_STRICT_THRESHOLD)
            allele_category[current_allele] = AlleleCategory::Adventageous;
        else if (advantageous_percent[current_allele] >= CATEGORY_INCLUSION_RELAXED_THRESHOLD)
            allele_category[current_allele] = AlleleCategory::Possible_advantageous;
        
        //unassigned_percent[current_allele] = unassigned_counter[current_allele] * 100 / _result_vector.size();
        
        // last in order to override any fitness value assignments
        if (current_allele == wt_allele) {
            allele_category[current_allele] = AlleleCategory::WT;
        }
    }
    
                                                                                                            
    // get settings for min and max fitness limits for prior
    /*
    MATRIX_TYPE fitness_limits_matrix(2, _num_alleles);
    for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
        std::string current_allele_fitness_str = fits_constants::PARAM_ALLELE_FITNESS + std::to_string(current_allele);
        std::string current_allele_min_fitness = fits_constants::PARAM_ALLELE_MIN_FITNESS + std::to_string(current_allele);
        std::string current_allele_max_fitness = fits_constants::PARAM_ALLELE_MAX_FITNESS + std::to_string(current_allele);
        
        fitness_limits_matrix(0, current_allele) = _zparams.GetFloat(current_allele_min_fitness, -1.0);
        fitness_limits_matrix(1, current_allele) = _zparams.GetFloat(current_allele_max_fitness, -1.0);
        
        if ( fitness_limits_matrix(0, current_allele) < 0  ||
             fitness_limits_matrix(0, current_allele) < 0 ) {
            
            throw "Missing fitness value for allele " + std::to_string(current_allele);
        }
    }
     */
    
    // auto num_prior_samples = _zparams.GetInt(fits_constants::PARAM_SIM_REPEATS);
    
    for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
        
        if ( current_allele == wt_allele ) {
            levenes_pval[current_allele] = -1.0f;
            continue;
        }
        // generate prior
        /*
        auto current_allele_min_fitness = fitness_limits_matrix(0, current_allele);
        auto current_allele_max_fitness = fitness_limits_matrix(1, current_allele);
        
        if ( current_allele_min_fitness == current_allele_max_fitness ) {
            std::cout << "skipping allele " << current_allele << std::endl;
            continue;
        }
        std::vector<FLOAT_TYPE> min_fitness_vec( num_prior_samples, current_allele_min_fitness );
        std::vector<FLOAT_TYPE> max_fitness_vec( num_prior_samples, current_allele_max_fitness );
        
        
        
        PriorDistributionType prior_type = PriorDistributionType::UNDEFINED;
        
        if ( _prior_distrib.compare(fits_constants::PARAM_PRIOR_DISTRIB_UNIFORM) == 0 ) {
            prior_type = PriorDistributionType::UNIFORM;
        }
        if ( _prior_distrib.compare(fits_constants::PARAM_PRIOR_DISTRIB_COMPOSITE) == 0 ) {
            prior_type = PriorDistributionType::FITNESS_COMPOSITE;
        }
        if ( _prior_distrib.compare(fits_constants::PARAM_PRIOR_DISTRIB_SMOOTHED_COMPOSITE) == 0 ) {
            prior_type = PriorDistributionType::SMOOTHED_COMPOSITE;
        }
        if ( _prior_distrib.compare(fits_constants::PARAM_PRIOR_DISTRIB_LOGNORMAL) == 0 ) {
            prior_type = PriorDistributionType::FITNESS_LOGNORMAL;
        }
        if ( prior_type == UNDEFINED ) {
            throw "ResultStats: prior not defined correctly";
        }
        
        
        PriorSampler<FLOAT_TYPE> sampler(min_fitness_vec, max_fitness_vec, PriorDistributionType::FITNESS_COMPOSITE);
        
        auto prior_vec = sampler.SamplePrior(1)[0];
        
        //std::cout << "prior size=" << prior_vec.size() << std::endl;
        */
        
        std::vector<FLOAT_TYPE> current_allele_prior_vec;
        for ( auto prior_vec : _prior_distrib ) {
            current_allele_prior_vec.push_back( prior_vec[current_allele] );
            // std::cout << "pushed to prior vec: " << prior_vec[current_allele] << std::endl;
        }
        
        boost::numeric::ublas::matrix_column<MATRIX_TYPE> current_allele_posterior_col( posterior_matrix, current_allele );
        
        std::vector<FLOAT_TYPE> current_allele_posterior_vec;
        current_allele_posterior_vec.resize(current_allele_posterior_col.size());
        
        std::copy( current_allele_posterior_col.cbegin(),
                   current_allele_posterior_col.cend(),
                  current_allele_posterior_vec.begin() );
        
        //auto partial_allele_prior_vec = whole_allele_prior_vec;
        //auto partial_allele_prior_vec = DownsampleVector( whole_allele_prior_vec, current_allele_posterior_vec.size() );
        //std::cout << "size of partial prior vec: " << current_allele_prior_vec.size() << "  posterior vec: " << current_allele_posterior_vec.size() << std::endl;
        FLOAT_TYPE levenes_p = -1.0f;
        try {
            levenes_p = LevenesTest2( current_allele_posterior_vec, current_allele_prior_vec );
        }
        catch ( const char* str ) {
            std::string tmp_str = "Error while calculating Levene's test: " + std::string(str);
            throw tmp_str;
        }
        catch ( std::string str ) {
            std::string tmp_str = "Error while calculating Levene's test: " + str;
            throw tmp_str;
        }
        catch (...) {
            std::string tmp_str = "Unknown error while calculating Levene's test.";
            throw tmp_str;
        }
        
        levenes_pval[current_allele] = levenes_p;
    }
}


std::string ResultsStats::GetSummaryFitness( bool table_only )
{
    
    std::stringstream ss;
    
    /*
    ss << "Quick Report - Fitness" << std::endl;
    ss << "=============================" << std::endl;

    ss << boost::format("%-10s") % "allele";
    ss << boost::format("%-10s") % "fitness";
    ss << boost::format("%-10s") % "category";
    ss << std::endl;

    for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {

     
         // Changed 2018-07-03
         // Adding insignificant pval to the warning
     
        if ( _allele_Nu[current_allele] >= 1.0 &&
            levenes_pval[current_allele] < fits_constants::LEVENES_SIGNIFICANCE ) {
            ss << boost::format(" %-10d") % current_allele;
        }
        else {
            ss << boost::format(" *%-9d") % current_allele;
            //Nu_flag = true;
        }

        ss << boost::format("%-10.3d") % allele_median_fitness[current_allele];
        ss << boost::format("%-10s") % AlleleCategory2String(allele_category[current_allele]);

        ss << std::endl;
    }

    ss << std::endl;
    */

    if ( !table_only ) {
        ss << "Fitness Report" << std::endl;
        ss << "===============" << std::endl;
        ss << GetSummaryHeader();
        
        
        auto tmp_size = _zparams.GetInt(fits_constants::PARAM_POPULATION_SIZE, -1);
        ss << "Population size (N) is " << tmp_size;
        
        if ( _zparams.GetInt(fits_constants::PARAM_SAMPLE_SIZE, 0) > 0 ) {
            ss << " (sampled " << _zparams.GetInt(fits_constants::PARAM_SAMPLE_SIZE, 0) << ")";
        }
        ss << std::endl;
        
        
        auto tmp_scaling_str = _zparams.GetString( fits_constants::PARAM_SCALING,
                                                  fits_constants::PARAM_SCALING_DEFAULT );
        
        if ( tmp_scaling_str.compare(fits_constants::PARAM_SCALING_OFF) == 0 ) {
            // ss << "Data has not been scaled" << std::endl;
        }
        else {
            ss << "Data was scaled using " << tmp_scaling_str << std::endl;
        }
        
        
        if ( _single_mutrate_used ) {
            ss << "Used a single mutation rate." << std::endl;
        }
        else {
            ss << "Used individual mutation rates." << std::endl;
        }
        
        ss << "Distance metric: " << _distance_metric << std::endl;
        
    }
    
    
    // warning for Nu
    std::vector<int> bad_alleles_vec;
    for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
        if ( _allele_Nu[current_allele] < 1.0f ) {
            bad_alleles_vec.push_back(current_allele);
        }
    }
    
    if ( !table_only ) {
        auto prior_str = _zparams.GetString( fits_constants::PARAM_FITNESS_PRIOR_DISTRIB,
                                          fits_constants::PARAM_PRIOR_DISTRIB_FITNESS_DEFAULT );
        
        ss << "Prior used: " << prior_str << std::endl;
        
        if ( !bad_alleles_vec.empty() ) {
            ss << "** WARNING: inference may be unreliable (Nu<1) for alleles (*): ";
            
            for ( auto allele_num : bad_alleles_vec ) {
                ss << allele_num << " ";
            }
            
            ss << " **" << std::endl << std::endl;
        }
        
        if ( _rejection_threshold > 0.0f ) {
            ss << "Rejection threshold set to " << boost::format("%-10d") % _rejection_threshold << std::endl;
        }
        
        ss << "--------------------" << std::endl;
    }
    
    // allele_MAD
    
    // TABLE PRINTING
    
    ss << boost::format("%-10s") % "allele";
    ss << boost::format("%-10s") % "median";
    ss << boost::format("%-10s") % "MAD";
    // ss << boost::format("%-10s") % "mean";
    ss << boost::format("%-10s") % "low";
    ss << boost::format("%-10s") % "high";
    ss << boost::format("%-10s") % "DEL(%)";
    ss << boost::format("%-10s") % "NEU(%)";
    ss << boost::format("%-10s") % "ADV(%)";
    ss << boost::format("%-10s") % "category";
    //ss << boost::format("%-10s") % "minldist";
    //ss << boost::format("%-10s") % "maxldist";
    //
    ss << boost::format("%-10s") % "pval";
    ss << std::endl;
    
    //bool Nu_flag = false;
    for ( auto current_allele=0; current_allele<_num_alleles; ++current_allele ) {
        
        /*
         Changed 2018-07-03
         Adding insignificant pval to the warning
         */
        if ( _allele_Nu[current_allele] >= 1.0 &&
            levenes_pval[current_allele] < fits_constants::LEVENES_SIGNIFICANCE ) {
            ss << boost::format("%-10d") % current_allele;
            //ss << boost::format("%.2e") % current_allele;
            
        }
        else {
            ss << boost::format("*%-9d") % current_allele;
            //ss << boost::format("*%.2e") % current_allele;
            
            //Nu_flag = true;
        }
        
        ss << boost::format("%-10.3f") % allele_median_fitness[current_allele];
        ss << boost::format("%-10.3f") % allele_MAD[current_allele];
         // ss << boost::format("%-10.3f") % allele_mean_fitness[current_allele];
        ss << boost::format("%-10.3f") % allele_min_fitness[current_allele];
        ss << boost::format("%-10.3f") % allele_max_fitness[current_allele];
        ss << boost::format("%-10.1f") % deleterious_percent[current_allele];
        ss << boost::format("%-10.1f") % neutral_percent[current_allele];
        ss << boost::format("%-10.1f") % advantageous_percent[current_allele];
        ss << boost::format("%-10s") % AlleleCategory2String(allele_category[current_allele]);
        //ss << boost::format("%-10.3d") % _distance_min;
        //ss << boost::format("%-10.3d") % _distance_max;
        //ss << boost::format("%-10.3d") % allele_MAD[current_allele];
        
        if ( allele_category[current_allele] == WT ) {
            ss << boost::format("%-10s") % "N/A";
        }
        else {
            //ss << boost::format("%-10.3d") % levenes_pval[current_allele];
            ss << boost::format("%.2e") % levenes_pval[current_allele];
            
        }
        
        ss << std::endl;
    }
    
    if ( !table_only ) {
        if ( _zparams.GetInt( fits_constants::PARAM_DUMP_PARAMETERS, 0) > 1 ) {
            ss << _zparams.GetAllParameters() << std::endl;
        }
    }
    
    return ss.str();
}


std::string ResultsStats::AlleleCategory2String( AlleleCategory category )
{
    std::string tmp_str = "";
    
    switch ( category ) {
            
        case AlleleCategory::Undefined:
            tmp_str = "???";
            break;
            
        case AlleleCategory::Adventageous:
            tmp_str = "ADV";
            break;
            
        case AlleleCategory::Possible_advantageous:
            tmp_str = "?ADV";
            break;
            
        case AlleleCategory::Deleterious:
            tmp_str = "DEL";
            break;
            
        case AlleleCategory::Possible_deleterious:
            tmp_str = "?DEL";
            break;
            
        case AlleleCategory::Neutral:
            tmp_str = "NEU";
            break;
        case AlleleCategory::Possible_neutral:
            tmp_str = "?NEU";
            break;
            
        case AlleleCategory::Lethal:
            tmp_str = "LTH";
            break;
            
        case AlleleCategory::Possible_lethal:
            tmp_str = "?LTH";
            break;
            
        case AlleleCategory::WT:
            tmp_str = "WT";
            break;
    }
    
    return tmp_str;
}

