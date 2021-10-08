
#include "GenericReconEigenChannelGadget.h"
#include <iomanip>

#include "hoNDArray_reductions.h"
#include "mri_core_def.h"

namespace Gadgetron {

    GenericReconEigenChannelGadget::GenericReconEigenChannelGadget() : BaseClass()
    {
    }

    GenericReconEigenChannelGadget::~GenericReconEigenChannelGadget()
    {
    }

    int GenericReconEigenChannelGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        try
        {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...)
        {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        if (!h.acquisitionSystemInformation)
        {
            GDEBUG("acquisitionSystemInformation not found in header. Bailing out");
            return GADGET_FAIL;
        }

        // -------------------------------------------------

        size_t NE = h.encoding.size();
        num_encoding_spaces_ = NE;
        GDEBUG_CONDITION_STREAM(verbose.value(), "Number of encoding spaces: " << NE);

        calib_mode_.resize(NE, ISMRMRD_noacceleration);

        data_KLT_.resize(NE);
        sms_ref_KLT_.resize(NE);
        ref_KLT_.resize(NE);

        for (size_t e = 0; e < h.encoding.size(); e++)
        {
            ISMRMRD::EncodingSpace e_space = h.encoding[e].encodedSpace;
            ISMRMRD::EncodingSpace r_space = h.encoding[e].reconSpace;
            ISMRMRD::EncodingLimits e_limits = h.encoding[e].encodingLimits;

            if (!h.encoding[e].parallelImaging)
            {
                GDEBUG_STREAM("Parallel Imaging section not found in header");
                calib_mode_[e] = ISMRMRD_noacceleration;
            }
            else
            {

                ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;
                std::string calib = *p_imaging.calibrationMode;

                bool separate = (calib.compare("separate") == 0);
                bool embedded = (calib.compare("embedded") == 0);
                bool external = (calib.compare("external") == 0);
                bool interleaved = (calib.compare("interleaved") == 0);
                bool other = (calib.compare("other") == 0);

                calib_mode_[e] = Gadgetron::ISMRMRD_noacceleration;
                if (p_imaging.accelerationFactor.kspace_encoding_step_1 > 1 || p_imaging.accelerationFactor.kspace_encoding_step_2 > 1)
                {
                    if (interleaved)
                        calib_mode_[e] = Gadgetron::ISMRMRD_interleaved;
                    else if (embedded)
                        calib_mode_[e] = Gadgetron::ISMRMRD_embedded;
                    else if (separate)
                        calib_mode_[e] = Gadgetron::ISMRMRD_separate;
                    else if (external)
                        calib_mode_[e] = Gadgetron::ISMRMRD_external;
                    else if (other)
                        calib_mode_[e] = Gadgetron::ISMRMRD_other;
                }
            }
        }

        return GADGET_OK;
    }

    int GenericReconEigenChannelGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1)
    {
        if (perform_timing.value()) { gt_timer_.start("GenericReconEigenChannelGadget::process"); }

        process_called_times_++;

        IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->rbit_.size() > num_encoding_spaces_)
        {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->rbit_.size() << " instead of " << num_encoding_spaces_);
        }

        // for every encoding space, prepare the recon_bit_->rbit_[e].ref_
        size_t e, n, s, slc;
        for (e = 0; e < recon_bit_->rbit_.size(); e++)
        {
            auto & rbit = recon_bit_->rbit_[e];
            std::stringstream os;
            os << "_encoding_" << e;

            hoNDArray< std::complex<float> >& data = recon_bit_->rbit_[e].data_.data_;

            size_t RO = data.get_size(0);
            size_t E1 = data.get_size(1);
            size_t E2 = data.get_size(2);
            size_t CHA = data.get_size(3);
            size_t N = data.get_size(4);
            size_t S = data.get_size(5);
            size_t SLC = data.get_size(6);

            GDEBUG_STREAM("GenericReconEigenChannelGadget - incoming data array : [RO E1 E2 CHA N S SLC] - [" << RO << " " << E1 << " " << E2 << " " << CHA << " " << N << " " << S << " " << SLC << "]");

            if(data.get_number_of_elements()==0)
            {
                m1->release();
                return GADGET_OK;
            }

            // whether it is needed to update coefficients
            bool recompute_coeff_data = false;
            if ( (data_KLT_[e].size()!=SLC) || update_eigen_channel_coefficients.value() )
            {
                recompute_coeff_data = true;
            }
            else
            {
                if(data_KLT_[e].size() == SLC)
                {
                    for (slc = 0; slc < SLC; slc++)
                    {
                        if (data_KLT_[e][slc].size() != S)
                        {
                            recompute_coeff_data = true;
                            break;
                        }
                        else
                        {
                            for (s = 0; s < S; s++)
                            {
                                if (data_KLT_[e][slc][s].size() != N)
                                {
                                    recompute_coeff_data = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            bool recompute_coeff_sms_ref = false;
            if ( (sms_ref_KLT_[e].size()!=SLC) || update_eigen_channel_coefficients.value() )
            {
                recompute_coeff_sms_ref = true;
            }
            else
            {
                if(sms_ref_KLT_[e].size() == SLC)
                {
                    for (slc = 0; slc < SLC; slc++)
                    {
                        if (sms_ref_KLT_[e][slc].size() != S)
                        {
                            recompute_coeff_sms_ref = true;
                            break;
                        }
                        else
                        {
                            for (s = 0; s < S; s++)
                            {
                                if (sms_ref_KLT_[e][slc][s].size() != N)
                                {
                                    recompute_coeff_sms_ref = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            bool recompute_coeff_ref = false;
            if ( (sms_ref_KLT_[e].size()!=SLC) || update_eigen_channel_coefficients.value() )
            {
                recompute_coeff_ref = true;
            }
            else
            {
                if(sms_ref_KLT_[e].size() == SLC)
                {
                    for (slc = 0; slc < SLC; slc++)
                    {
                        if (sms_ref_KLT_[e][slc].size() != S)
                        {
                            recompute_coeff_ref = true;
                            break;
                        }
                        else
                        {
                            for (s = 0; s < S; s++)
                            {
                                if (sms_ref_KLT_[e][slc][s].size() != N)
                                {
                                    recompute_coeff_ref = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            bool average_N = average_all_ref_N.value();
            bool average_S = average_all_ref_S.value();


            if (rbit.ref_) {
                if (recompute_coeff_ref) {
                    Gadgetron::compute_eigen_channel_coefficients(
                        rbit.ref_->data_, average_N, average_S, (calib_mode_[e] == Gadgetron::ISMRMRD_interleaved), N,
                        S, upstream_coil_compression_thres.value(), upstream_coil_compression_num_modesKept.value(),
                        ref_KLT_[e]);
                }

                if (verbose.value())
                {
                    hoNDArray< std::complex<float> > E;

                    for (slc = 0; slc < SLC; slc++)
                    {
                        for (s = 0; s < S; s++)
                        {
                            for (n = 0; n < N; n++)
                            {
                                ref_KLT_[e][slc][s][n].eigen_value(E);
                            }
                        }
                    }
                } else {
                    if(average_N && average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, SLC : " << slc << " - " << ref_KLT_[e][slc][0][0].output_length() << " out of " << CHA);
                        }
                    }
                    else if(average_N && !average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (s = 0; s < S; s++)
                            {
                                GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, [SLC S] : [" << slc << " " << s << "] - " << ref_KLT_[e][slc][s][0].output_length() << " out of " << CHA);
                            }
                        }
                    }
                    else if(!average_N && average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (n = 0; n < N; n++)
                            {
                                GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, [SLC N] : [" << slc << " " << n << "] - " << ref_KLT_[e][slc][0][n].output_length() << " out of " << CHA);
                            }
                        }
                    }
                    else if(!average_N && !average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (s = 0; s < S; s++)
                            {
                                for (n = 0; n < N; n++)
                                {
                                    GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, [SLC S N] : [" << slc << " " << s << " " << n << "] - " << ref_KLT_[e][slc][s][n].output_length() << " out of " << CHA);
                                }
                            }
                        }
                    }
                }


                if (!debug_folder_full_path_.empty()){
                    gt_exporter_.export_array_complex(rbit.ref_->data_, debug_folder_full_path_ + "ref_before_KLT" + os.str());
                }

                Gadgetron::apply_eigen_channel_coefficients(ref_KLT_[e], rbit.ref_->data_);

                if (!debug_folder_full_path_.empty()){
                    gt_exporter_.export_array_complex(rbit.ref_->data_, debug_folder_full_path_ + "ref_after_KLT" + os.str());
                }
            }

            if (rbit.data_.data_.get_number_of_elements() > 0) {
                if (recompute_coeff_data) {
                    Gadgetron::compute_eigen_channel_coefficients(
                        rbit.data_.data_, average_N, average_S, (calib_mode_[e] == Gadgetron::ISMRMRD_interleaved), N,
                        S, upstream_coil_compression_thres.value(), upstream_coil_compression_num_modesKept.value(),
                        data_KLT_[e]);
                }

                if (verbose.value()) {
                    hoNDArray<std::complex<float>> E;

                    for (slc = 0; slc < SLC; slc++) {
                        for (s = 0; s < S; s++) {
                            for (n = 0; n < N; n++) {
                                data_KLT_[e][slc][s][n].eigen_value(E);
                            }
                        }
                    }
                } else {
                    if(average_N && average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, SLC : " << slc << " - " << data_KLT_[e][slc][0][0].output_length() << " out of " << CHA);
                        }
                    }
                    else if(average_N && !average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (s = 0; s < S; s++)
                            {
                                GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, [SLC S] : [" << slc << " " << s << "] - " << data_KLT_[e][slc][s][0].output_length() << " out of " << CHA);
                            }
                        }
                    }
                    else if(!average_N && average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (n = 0; n < N; n++)
                            {
                                GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, [SLC N] : [" << slc << " " << n << "] - " << data_KLT_[e][slc][0][n].output_length() << " out of " << CHA);
                            }
                        }
                    }
                    else if(!average_N && !average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (s = 0; s < S; s++)
                            {
                                for (n = 0; n < N; n++)
                                {
                                    GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, [SLC S N] : [" << slc << " " << s << " " << n << "] - " << data_KLT_[e][slc][s][n].output_length() << " out of " << CHA);
                                }
                            }
                        }
                    }
                }

                if (!debug_folder_full_path_.empty()){
                    gt_exporter_.export_array_complex(rbit.data_.data_, debug_folder_full_path_ + "data_before_KLT" + os.str());
                }

                // apply KL coefficients
                Gadgetron::apply_eigen_channel_coefficients(data_KLT_[e], rbit.data_.data_);

                if (!debug_folder_full_path_.empty()){
                    gt_exporter_.export_array_complex(rbit.data_.data_, debug_folder_full_path_ + "data_after_KLT" + os.str());
                }
            }

            if (rbit.sms_ref_) {
                if (recompute_coeff_sms_ref) {
                    Gadgetron::compute_eigen_channel_coefficients(
                        rbit.sms_ref_->data_, average_N, average_S, (calib_mode_[e] == Gadgetron::ISMRMRD_interleaved),
                        N, S, upstream_coil_compression_thres.value(), upstream_coil_compression_num_modesKept.value(),
                        sms_ref_KLT_[e]);
                }

                if (verbose.value())
                {
                    hoNDArray< std::complex<float> > E;

                    for (slc = 0; slc < SLC; slc++)
                    {
                        for (s = 0; s < S; s++)
                        {
                            for (n = 0; n < N; n++)
                            {
                                sms_ref_KLT_[e][slc][s][n].eigen_value(E);
                            }
                        }
                    }
                } else {
                    if(average_N && average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, SLC : " << slc << " - " << sms_ref_KLT_[e][slc][0][0].output_length() << " out of " << CHA);
                        }
                    }
                    else if(average_N && !average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (s = 0; s < S; s++)
                            {
                                GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, [SLC S] : [" << slc << " " << s << "] - " << sms_ref_KLT_[e][slc][s][0].output_length() << " out of " << CHA);
                            }
                        }
                    }
                    else if(!average_N && average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (n = 0; n < N; n++)
                            {
                                GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, [SLC N] : [" << slc << " " << n << "] - " << sms_ref_KLT_[e][slc][0][n].output_length() << " out of " << CHA);
                            }
                        }
                    }
                    else if(!average_N && !average_S)
                    {
                        for (slc = 0; slc < SLC; slc++)
                        {
                            for (s = 0; s < S; s++)
                            {
                                for (n = 0; n < N; n++)
                                {
                                    GDEBUG_STREAM("GenericReconEigenChannelGadget - Number of modes kept, [SLC S N] : [" << slc << " " << s << " " << n << "] - " << sms_ref_KLT_[e][slc][s][n].output_length() << " out of " << CHA);
                                }
                            }
                        }
                    }
                }

                if (!debug_folder_full_path_.empty()){
                    gt_exporter_.export_array_complex(rbit.sms_ref_->data_, debug_folder_full_path_ + "sms_ref_before_KLT" + os.str());
                }

                Gadgetron::apply_eigen_channel_coefficients(sms_ref_KLT_[e], rbit.sms_ref_->data_);

                if (!debug_folder_full_path_.empty()){
                    gt_exporter_.export_array_complex(rbit.sms_ref_->data_, debug_folder_full_path_ + "sms_ref_after_KLT" + os.str());
                }
            }
        }

        if (perform_timing.value()) { gt_timer_.stop(); }

        if (this->next()->putq(m1) < 0)
        {
            GERROR_STREAM("Put IsmrmrdReconData to Q failed ... ");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(GenericReconEigenChannelGadget)
}
