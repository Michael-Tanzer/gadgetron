/** \file   SMSReconCartesianGrappaGadget.h
    \brief
    \author Michael Tanzer
*/

#pragma once

#include "GenericReconGadget.h"


namespace Gadgetron {

    class EXPORTGADGETSMRICORE SMSReconCartesianGrappaGadget : public GenericReconGadget
    {
    public:
        GADGET_DECLARE(SMSReconCartesianGrappaGadget);

        typedef GenericReconGadget BaseClass;
//        typedef Gadgetron::GenericReconCartesianGrappaObj< std::complex<float> > ReconObjType;

        SMSReconCartesianGrappaGadget();
        ~SMSReconCartesianGrappaGadget() override;

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------
    protected:

        // --------------------------------------------------
        // variable for recon
        // --------------------------------------------------
        // record the recon kernel, coil maps etc. for every encoding space
//        std::vector< ReconObjType > recon_obj_;

        // --------------------------------------------------
        // gadget functions
        // --------------------------------------------------
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);
 };
}
