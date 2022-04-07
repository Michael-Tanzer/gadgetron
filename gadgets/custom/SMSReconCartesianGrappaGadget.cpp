
#include "SMSReconCartesianGrappaGadget.h"
#include "custom_mri_core_grappa.h"
#include "mri_core_grappa.h"
#include "hoNDArray_reductions.h"
#include <any>
#include "GenericReconCartesianGrappaGadget.h"


namespace Gadgetron {

    SMSReconCartesianGrappaGadget:: SMSReconCartesianGrappaGadget() : BaseClass() {
    }

    SMSReconCartesianGrappaGadget::~ SMSReconCartesianGrappaGadget() {
    }

    int SMSReconCartesianGrappaGadget::process_config(ACE_Message_Block *mb) {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);
        GDEBUG("DEBUGGING CUSTOM SMS GRAPPA  ----------------------------  1\n");

        // -------------------------------------------------

        ISMRMRD::IsmrmrdHeader h;
        try {
            deserialize(mb->rd_ptr(), h);
        }
        catch (...) {
            GDEBUG("Error parsing ISMRMRD Header");
        }

        return GADGET_OK;
    }

//    int SMSReconCartesianGrappaGadget::process(Gadgetron::GadgetContainerMessage<IsmrmrdReconData> *m1) {
//        GDEBUG("DEBUGGING SMS GRAPPA  ----------------------------  2\n");
//    }
    int SMSReconCartesianGrappaGadget::process(Gadgetron::GadgetContainerMessage<IsmrmrdReconData> *m1) {
        process_called_times_++;

        GDEBUG("DEBUGGING SMS GRAPPA  ----------------------------  2\n");

        IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
        if (recon_bit_->rbit_.size() > num_encoding_spaces_) {
            GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : "
                         << recon_bit_->rbit_.size() << " instead of " << num_encoding_spaces_);
        }

        GadgetContainerMessage<std::vector<ISMRMRD::Waveform>>* wav =
            AsContainerMessage<std::vector<ISMRMRD::Waveform>>(m1->cont());
        if (wav) {
            if (verbose.value()) {
                GDEBUG_STREAM("Incoming recon_bit with " << wav->getObjectPtr()->size() << " wave form samples ");
            }
        }

        // for every encoding space
        for (size_t e = 0; e < recon_bit_->rbit_.size(); e++) {
            std::stringstream os;
            os << "_encoding_" << e << "_" << process_called_times_;

            GDEBUG_CONDITION_STREAM(verbose.value(),
                                    "Calling " << process_called_times_ << " , encoding space : " << e);
            GDEBUG_CONDITION_STREAM(verbose.value(),
                                    "======================================================================");

            // ---------------------------------------------------------------
            // export incoming data
            GDEBUG("DEBUGGING CUSTOM GRAPPA  ----------------------------  3\n");
            GDEBUG_STREAM("FILEPATH" << debug_folder_full_path_ << "\n");
            if (!debug_folder_full_path_.empty()) {
                gt_exporter_.export_array_complex(recon_bit_->rbit_[e].data_.data_,
                                                  debug_folder_full_path_ + "data" + os.str());
            }

            if (!debug_folder_full_path_.empty() && recon_bit_->rbit_[e].data_.trajectory_) {
                if (recon_bit_->rbit_[e].ref_->trajectory_->get_number_of_elements() > 0) {
                    gt_exporter_.export_array(*(recon_bit_->rbit_[e].data_.trajectory_),
                                              debug_folder_full_path_ + "data_traj" + os.str());
                }
            }

            // ---------------------------------------------------------------
            auto casted_recon_obj_ = std::any_cast<Gadgetron::GenericReconCartesianGrappaObj<std::complex<float>>>(recon_bit_->rbit_[e].additional_data.at("grappa_recon_obj"));

            if (recon_bit_->rbit_[e].sms_ref_) {
            }

//            try {
//                auto casted_recon_obj_ = std::any_cast<GenericReconCartesianGrappaObj<std::complex<float>>>(recon_bit_->rbit_[e].additional_data.at("grappa_recon_obj"));
//            }
//            catch (const std::out_of_range&) {
//                GERROR("The SMS gadget needs to be places after the GenericReconCartesiaGrappaGadget");
//            }

            GDEBUG("HELLO");


//            if (recon_bit_->rbit_[e].ref_) {
//
//                // this was the very last step
//                recon_bit_->rbit_[e].ref_ = Core::none;
//            }

//            if (recon_bit_->rbit_[e].data_.data_.get_number_of_elements() > 0) {
//
//                // pass down waveform  - WHAT DOES THIS DO?
//                if (wav)
//                    recon_obj_[e].recon_res_.waveform_ = *wav->getObjectPtr();
//                recon_obj_[e].recon_res_.acq_headers_ = recon_bit_->rbit_[e].data_.headers_;
//
//                // NORMAL GRAPPA SENDS OUT IMAGES? IS THIS RIGHT? WHAT DOES THIS GADGET GET AS INPUT?
//                this->send_out_image_array(recon_obj_[e].recon_res_, e, image_series.value() + ((int)e + 1),
//                                           GADGETRON_IMAGE_REGULAR);
//            }
//
//            // this is clearing everything. I suspect I need to change the grappa gadget so that it sends all the data and not images.
//            recon_obj_[e].recon_res_.data_.clear();
//            recon_obj_[e].gfactor_.clear();
//            recon_obj_[e].recon_res_.headers_.clear();
//            recon_obj_[e].recon_res_.meta_.clear();
        }

        // needs to be here apparently
        m1->release();

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(SMSReconCartesianGrappaGadget)
}