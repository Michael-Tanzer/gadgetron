#include "CustomEPIReconXGadget.h"
#include "ismrmrd/xml.h"

#ifdef USE_OMP
#include "omp.h"
#endif // USE_OMP

namespace Gadgetron{

  CustomEPIReconXGadget::CustomEPIReconXGadget() {}
  CustomEPIReconXGadget::~CustomEPIReconXGadget() {}

int CustomEPIReconXGadget::process_config(ACE_Message_Block* mb)
{
  ISMRMRD::IsmrmrdHeader h;
  ISMRMRD::deserialize(mb->rd_ptr(),h);
  
  
  verboseMode_ = verboseMode.value();
  // todo: add other attributes here

  if (h.encoding.size() == 0) {
    GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
    GDEBUG("This Gadget needs an encoding description\n");
    return GADGET_FAIL;
  }

  GDEBUG("Number of encoding spaces = %d\n", h.encoding.size());

  // Get the encoding space and trajectory description
  ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
  ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
  ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;
  ISMRMRD::TrajectoryDescription traj_desc;

  if (h.encoding[0].trajectoryDescription) {
    traj_desc = *h.encoding[0].trajectoryDescription;
  } else {
    GDEBUG("Trajectory description missing");
    return GADGET_FAIL;
  }

  if (traj_desc.identifier != "ConventionalEPI") {
    GDEBUG("Expected trajectory description identifier 'ConventionalEPI', not found.");
    return GADGET_FAIL;
  }

  // Primary encoding space is for EPI
  reconx.encodeNx_  = e_space.matrixSize.x;
  reconx.encodeFOV_ = e_space.fieldOfView_mm.x;
  reconx.reconNx_   = r_space.matrixSize.x;
  reconx.reconFOV_  = r_space.fieldOfView_mm.x;

  GDEBUG("1 - ENCODE NX: %d\n", reconx.encodeNx_);
  GDEBUG("1 - ENCODE FOV: %d\n", reconx.encodeFOV_);
  GDEBUG("1 - RECON NX: %d\n", reconx.reconNx_);
  GDEBUG("1 - RECON FOV: %d\n", reconx.reconFOV_);

  if (reconx.reconNx_ == 0) {
    GDEBUG("IN IF NUMBER 1\n");
    reconx.reconNx_ = 160;
  }

  GDEBUG("2 - RECON NX: %d\n", reconx.reconNx_);
  
  // TODO: we need a flag that says it's a balanced readout.
  for (std::vector<ISMRMRD::UserParameterLong>::iterator i (traj_desc.userParameterLong.begin()); i != traj_desc.userParameterLong.end(); ++i) {
    if (i->name == "rampUpTime") {
      GDEBUG("IN IF NUMBER 2\n");
      reconx.rampUpTime_ = i->value;
    } else if (i->name == "rampDownTime") {
      GDEBUG("IN IF NUMBER 3\n");
      reconx.rampDownTime_ = i->value;
    } else if (i->name == "flatTopTime") {
      GDEBUG("IN IF NUMBER 4\n");
      reconx.flatTopTime_ = i->value;
    } else if (i->name == "acqDelayTime") {
      GDEBUG("IN IF NUMBER 5\n");
      reconx.acqDelayTime_ = i->value;
    } else if (i->name == "numSamples") {
      GDEBUG("IN IF NUMBER 6\n");
      reconx.numSamples_ = i->value;
      //reconx.numSamples_ = 160;
    }
  }
  GDEBUG("3 - RAMPUPTIME: %d\n", reconx.rampUpTime_);
  GDEBUG("3 - RAMPDOWNTIME: %d\n", reconx.rampDownTime_);
  GDEBUG("3 - FLATTOPTIME: %d\n", reconx.flatTopTime_);
  GDEBUG("3 - ACQDELAYTIME: %d\n", reconx.acqDelayTime_);
  GDEBUG("3 - NUMSAMPLES: %d\n", reconx.numSamples_);


  for (std::vector<ISMRMRD::UserParameterDouble>::iterator i (traj_desc.userParameterDouble.begin()); i != traj_desc.userParameterDouble.end(); ++i) {
    if (i->name == "dwellTime") {
      GDEBUG("IN IF NUMBER 7\n");
      reconx.dwellTime_ = i->value;
    }
  }
  GDEBUG("4 - DWELLTIME: %d\n", reconx.dwellTime_);


  // If the flat top time is not set in the header, then we assume that rampSampling is off
  // and we set the flat top time from the number of samples and the dwell time.
  if (reconx.flatTopTime_ == 0) {
      GDEBUG("IN IF NUMBER 8\n");
      reconx.flatTopTime_ = reconx.dwellTime_ * reconx.numSamples_;
  }
  GDEBUG("5 - FLATTOPTIME: %d\n", reconx.flatTopTime_);

  // Compute the trajectory
  reconx.computeTrajectory();

  // Second encoding space is an even readout for PAT REF e.g. FLASH
  if ( h.encoding.size() > 1 ) {
    GDEBUG("IN IF NUMBER 9\n");
    ISMRMRD::EncodingSpace e_space2 = h.encoding[1].encodedSpace;
    ISMRMRD::EncodingSpace r_space2 = h.encoding[1].reconSpace;
    reconx_other.encodeNx_  = r_space2.matrixSize.x;
    reconx_other.encodeFOV_ = r_space2.fieldOfView_mm.x;
    reconx_other.reconNx_   = r_space2.matrixSize.x;
    reconx_other.reconFOV_  = r_space2.fieldOfView_mm.x;
    reconx_other.numSamples_ = e_space2.matrixSize.x;
    oversamplng_ratio2_ = (float)e_space2.matrixSize.x / r_space2.matrixSize.x;
    reconx_other.dwellTime_ = 1.0;

    if (reconx_other.reconNx_ == 0) {
        reconx_other.reconNx_ = 160;
    }

    reconx_other.computeTrajectory();
  }
  GDEBUG("6 - ENCODE NX: %d\n", reconx_other.encodeNx_);
  GDEBUG("6 - ENCODE FOV: %d\n", reconx_other.encodeFOV_);
  GDEBUG("6 - RECON NX: %d\n", reconx_other.reconNx_);
  GDEBUG("6 - RECON FOV: %d\n", reconx_other.reconFOV_);
  GDEBUG("6 - NUMSAMPLES: %d\n", reconx_other.numSamples_);
  GDEBUG("6 - DEWLLTIME: %d\n", reconx_other.dwellTime_);

  GDEBUG("7 - RECON NX: %d\n", reconx_other.reconNx_);

#ifdef USE_OMP
  omp_set_num_threads(1);
#endif // USE_OMP

  return 0;
}

int CustomEPIReconXGadget::process(
          GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
      GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
{

  ISMRMRD::AcquisitionHeader hdr_in = *(m1->getObjectPtr());
  ISMRMRD::AcquisitionHeader hdr_out;
  hoNDArray<std::complex<float> > data_out;

  data_out.create(reconx.reconNx_, m2->getObjectPtr()->get_size(1));

  // Switch the reconstruction based on the encoding space (e.g. for FLASH Calibration)
  if (hdr_in.encoding_space_ref == 0) {
    GDEBUG("IN IF NUMBER 10\n");
    reconx.apply(*m1->getObjectPtr(), *m2->getObjectPtr(), hdr_out, data_out);
  }
  else
  {
    GDEBUG("IN IF NUMBER 11\n");
    if(reconx_other.encodeNx_>m2->getObjectPtr()->get_size(0)/ oversamplng_ratio2_)
    {
        GDEBUG("IN IF NUMBER 12\n");
        reconx_other.encodeNx_ = (int)(m2->getObjectPtr()->get_size(0) / oversamplng_ratio2_);
        GDEBUG("8 - ENCODE NX: %d", reconx_other.encodeNx_);
        reconx_other.computeTrajectory();
    }

    if (reconx_other.reconNx_>m2->getObjectPtr()->get_size(0) / oversamplng_ratio2_)
    {
        GDEBUG("IN IF NUMBER 13\n");
        reconx_other.reconNx_ = (int)(m2->getObjectPtr()->get_size(0) / oversamplng_ratio2_);
        GDEBUG("9 - RECON NX: %d", reconx_other.reconNx_);
    }

    if(reconx_other.numSamples_>m2->getObjectPtr()->get_size(0))
    {
        GDEBUG("IN IF NUMBER 14\n");
        reconx_other.numSamples_ = m2->getObjectPtr()->get_size(0);
        GDEBUG("14 - NUMSAMPLES: %d", reconx_other.numSamples_);
        reconx_other.computeTrajectory();
    }

    reconx_other.apply(*m1->getObjectPtr(), *m2->getObjectPtr(), hdr_out, data_out);
  }

  // Replace the contents of m1 with the new header and the contentes of m2 with the new data
  *m1->getObjectPtr() = hdr_out;
  *m2->getObjectPtr() = data_out;

  // It is enough to put the first one, since they are linked
  if (this->next()->putq(m1) == -1) {
    m1->release();
    GERROR("CustomEPIReconXGadget::process, passing data on to next gadget");
    return -1;
  }

  return 0;
}

GADGET_FACTORY_DECLARE(CustomEPIReconXGadget)
}


