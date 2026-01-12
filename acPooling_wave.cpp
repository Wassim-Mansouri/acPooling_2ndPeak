// This file is part of the HörTech Open Master Hearing Aid (openMHA)
// Copyright © 2015 2016 2018 2019 2020 2021 HörTech gGmbH
// Copyright © 2022 Hörzentrum Oldenburg gGmbH
//
// openMHA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, version 3 of the License.
//
// openMHA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License, version 3 for more details.
//
// You should have received a copy of the GNU Affero General Public License, 
// version 3 along with openMHA.  If not, see <http://www.gnu.org/licenses/>.

/*
 * This plugin receives a waveform through the ac space at each
 * iteration. Depending on the averaging method, the last consecutive
 * n frames are averaged and the computed average is send back to the
 * ac space. Also, the maximum of the computed average is written to
 * the ac space.
 */

#include "acPooling_wave.h"

#define PATCH_VAR(var) patchbay.connect(&var.valuechanged, this, &acPooling_wave::update_cfg)
#define INSERT_PATCH(var) insert_member(var); PATCH_VAR(var)

acPooling_wave_config::acPooling_wave_config(MHA_AC::algo_comm_t & ac,
                                             const mhaconfig_t in_cfg,
                                             acPooling_wave *_pooling):
    ac(ac),
    raw_p_name(_pooling->p_name.data),
    p(ac, _pooling->pool_name.data.c_str(), _pooling->numsamples.data, 1, false),
    p_biased(ac, _pooling->p_biased_name.data.c_str(), _pooling->numsamples.data, 1, false),
    like_ratio(ac, _pooling->like_ratio_name.data.c_str(), 1, 1, false),
    p_max(ac, _pooling->max_pool_ind_name.data.c_str(), 1, 1, false),
    p_second_max(ac, "pool_second_max", 1, 1, false), //  Second peak AC*
    pooling_ind(0),
    pooling_option(_pooling->pooling_type.data.get_index()),
    pooling_size(std::max(2u, (unsigned int)(_pooling->pooling_wndlen.data * in_cfg.srate / (in_cfg.fragsize * 1000)))),
    up_thresh(_pooling->upper_threshold.data),
    low_thresh(_pooling->lower_threshold.data),
    neigh(_pooling->neighbourhood.data),
    alpha(_pooling->alpha.data),
    pool(_pooling->numsamples.data, std::max(2u, (unsigned int)(_pooling->pooling_wndlen.data * in_cfg.srate / (in_cfg.fragsize * 1000)))),
    prob_bias_func(_pooling->prob_bias.data)
{
    pool.assign(0);
    p_max.assign((_pooling->numsamples.data - 1) / 2);
    //  second peak variables 
    second_max_ind = 0;
    second_max_value = 0;
    min_peak_separation = 2; //to be changed later.
    // 
}

acPooling_wave_config::~acPooling_wave_config() {}

void acPooling_wave_config::insert()
{
    p.insert();
    p_biased.insert();
    p_max.insert();
    like_ratio.insert();
    //  Inserted second peak 
    p_second_max.insert();
    p_second_max(0,0) = second_max_ind;
    //
}

mha_wave_t *acPooling_wave_config::process(mha_wave_t *wave)
{
    const mha_wave_t raw_p = MHA_AC::get_var_waveform(ac, raw_p_name.c_str());

    mha_real_t max = 0, sample_max, mean_p = 0;
    int max_ind = -1;

    for (unsigned int i = 0;  i < p.num_frames; ++i) {
        switch(pooling_option) {
        case 0:
            if (alpha > 0 )
                p(i, 0) -= pow(1 - alpha, pooling_size - 1) * alpha * pool.value(i, pooling_ind);
            else
                p(i, 0) -= pool.value(i, pooling_ind);
            break;
        case 2:
            if (alpha > 0)
                p(i, 0) -= pow(1 - alpha, pooling_size - 1) * alpha * pool.value(i, pooling_ind) / pooling_size;
            else
                p(i, 0) -= ( pool.value(i, pooling_ind) / pooling_size );
            break;
        }
    }

    for( unsigned int i = 0; i < p.num_frames; ++i ) {
        pool(i, pooling_ind) =  value(raw_p, i, 0);
        switch (pooling_option) {
        case 0:
            if (alpha > 0) 
                p(i, 0) = (1 - alpha) * p(i, 0) + alpha * pool.value(i, pooling_ind);
            else
                p(i, 0) += pool.value(i, pooling_ind);
            break;
        case 1:
            sample_max = 0;
            for (unsigned int j = 0; j < pooling_size; j++) {
                if (sample_max < pool.value(i, j))
                    sample_max = pool.value(i, j);
                p(i, 0) = sample_max;
            }
            break;
        case 2:
            if (alpha > 0)
                p(i, 0) = (1 - alpha) * p(i, 0) + alpha * pool.value(i, pooling_ind) / pooling_size;
            else
                p(i, 0) += ( pool.value(i, pooling_ind) / pooling_size );
            break;
        }

        mean_p += p.value(i, 0);
        p_biased(i, 0) = p.value(i, 0)*prob_bias_func(i, 0);

        // compute the max peaks with that seperation
        if (i == 0) {
            max = p_biased(i,0);
            max_ind = i;
            second_max_value = 0;
            second_max_ind = -1;
        } else {
            mha_real_t val = p_biased(i,0);
            if (val > max) {
                if (max_ind != -1 && std::abs((int)i - (int)max_ind) >= (int)min_peak_separation) {
                    second_max_value = max;
                    second_max_ind = max_ind;
                }
                max = val;
                max_ind = i;
            } else if (val > second_max_value) {
                if (std::abs((int)i - (int)max_ind) >= (int)min_peak_separation) {
                    second_max_value = val;
                    second_max_ind = i;
                }
            }
        }
        // 
    }

    if (max_ind == -1)
        max_ind = (int)((p.num_frames - 1) / 2);

    pooling_ind = (pooling_ind + 1)  % pooling_size;
    like_ratio(0, 0) = p_biased(max_ind, 0) / mean_p;

    if ( like_ratio(0, 0) >= up_thresh )
        p_max(0, 0) = max_ind;
    else if (like_ratio(0, 0) >= low_thresh) {
        if (neigh == -1)
            p_max(0, 0) = max_ind;
        else if (std::abs(max_ind - (int)p_max(0, 0)) <= neigh)
            p_max(0, 0) = max_ind;
    }

    insert();
    return wave;
}

acPooling_wave::acPooling_wave(MHA_AC::algo_comm_t & iac,
                               const std::string & configured_name)
    : MHAPlugin::plugin_t<acPooling_wave_config>("Pooling of several consecutive time frames",iac)
    , numsamples("This parameter determines the length of the wave to be pooled in samples", "37", "]0,]")
    , pooling_wndlen("This parameter determines the length of the pooling window in msec.", "300", "]0,]")
    , pooling_type("This parameter determines the pooling method applied to the pooling window.", "mean", "[max sum mean]")
    , upper_threshold("This parameter sets a threshold for finding the maximum probability. If the maximum is above this threshold, it is taken, even if it is not in the neighbourhood of the last estimated direction.", "0.75", "[0, 1]")
    , lower_threshold("This parameter sets a threshold for finding the minimum probability. If the maximum probability is below this threshold, the estimated direction of the last iteration is taken.", "0", "[0, 1]")
    , neighbourhood("This parameter defines the neighbourhood of the allowed change of the estimated direction between iterations. -1 means no neighbourhood.", "2", "[-1,]")
    , alpha("This parameter simulates the forgetting effect by weighting the frames within the pooling window,  e.g. p(n + 1) = (1 - alpha) * p(n) + alpha * p_new. 0 means no weighting.", "0.1", "[0, 1]")
    , p_name("The name of the AC variable of the frame, which is going to be pooled.", "p")
    , p_biased_name("The name of the AC variable of the biased frame, after pooling.", "prob_biased")
    , pool_name("The name of the AC variable of the averaged (pooled) frame.", "pool")
    , max_pool_ind_name("The name of the AC variable for the index of the maximum of the averaged frames", "pool_max"),
    second_max_pool_ind_name("The name of the AC variable for the index of the second peak", "pool_second_max")
    
    , like_ratio_name("The name of the AC variable for the likelihood ratios of the averaged frames", "like_ratio")
    , prob_bias("A multiplicative probability bias", "[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]")
{
    set_node_id(configured_name);

    INSERT_PATCH(numsamples);
    INSERT_PATCH(pooling_wndlen);
    INSERT_PATCH(pooling_type);
    INSERT_PATCH(upper_threshold);
    INSERT_PATCH(lower_threshold);
    INSERT_PATCH(neighbourhood);
    INSERT_PATCH(alpha);
    INSERT_PATCH(p_name);
    INSERT_PATCH(p_biased_name);
    INSERT_PATCH(pool_name);
    INSERT_PATCH(max_pool_ind_name);
    INSERT_PATCH(second_max_pool_ind_name);
    INSERT_PATCH(like_ratio_name);
    INSERT_PATCH(prob_bias);
}

acPooling_wave::~acPooling_wave() {}

void acPooling_wave::prepare(mhaconfig_t & signal_info)
{
    if( signal_info.domain != MHA_WAVEFORM )
        throw MHA_Error(__FILE__, __LINE__,
                        "This plug-in requires time-domain signals.");

    if (pooling_wndlen.data * signal_info.srate / 1000 < signal_info.fragsize * 2)
        throw MHA_Error (__FILE__, __LINE__,
                         "Pooling window size in samples (%f) is not allowed to be smaller "
                         "than the twice of the fragsize (%u).",
                         pooling_wndlen.data * signal_info.srate / 1000, signal_info.fragsize * 2);
    
    if (prob_bias.data.size() != static_cast<unsigned int>(numsamples.data))
        throw MHA_Error(__FILE__, __LINE__,
                        "prob_bias should have \"numsamples\" (%i) elements but has %zu elements",
                        numsamples.data, prob_bias.data.size());

    update_cfg();
    poll_config()->insert();
}

void acPooling_wave::update_cfg()
{
    if ( is_prepared() ) {
        acPooling_wave_config *config;
        config = new acPooling_wave_config( ac, input_cfg(), this );
        push_config( config );
    }
}

mha_wave_t * acPooling_wave::process(mha_wave_t * signal)
{
    return poll_config()->process( signal );
}

MHAPLUGIN_CALLBACKS(acPooling_wave,acPooling_wave,wave,wave)
MHAPLUGIN_DOCUMENTATION\
(acPooling_wave,
 "data-flow feature-extraction algorithm-communication",
        "This plugin computes an average over several consecutive frames using several different approaches (max, sum, mean, ...). Subsequently, the maximum of the average frame is delivered as well. The plugin receives the frames through the AC space and deliveres the averaged frame as well as its maximum to the AC space.\n"
        "...")

