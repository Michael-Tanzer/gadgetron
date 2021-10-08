
/** \file   mri_core_grappa.cpp
    \brief  GRAPPA implementation for 2D and 3D MRI parallel imaging
    \author Hui Xue

    References to the implementation can be found in:

    Griswold MA, Jakob PM, Heidemann RM, Nittka M, Jellus V, Wang J, Kiefer B, Haase A. 
    Generalized autocalibrating partially parallel acquisitions (GRAPPA). 
    Magnetic Resonance in Medicine 2002;47(6):1202-1210.

    Kellman P, Epstein FH, McVeigh ER. 
    Adaptive sensitivity encoding incorporating temporal filtering (TSENSE). 
    Magnetic Resonance in Medicine 2001;45(5):846-852.

    Breuer FA, Kellman P, Griswold MA, Jakob PM. .
    Dynamic autocalibrated parallel imaging using temporal GRAPPA (TGRAPPA). 
    Magnetic Resonance in Medicine 2005;53(4):981-985.

    Saybasili H., Kellman P., Griswold MA., Derbyshire JA. Guttman, MA. 
    HTGRAPPA: Real-time B1-weighted image domain TGRAPPA reconstruction. 
    Magnetic Resonance in Medicine 2009;61(6): 1425-1433. 
*/

#include "mri_core_grappa.h"
#include "mri_core_utility.h"
#include "hoMatrix.h"
#include "hoNDArray_linalg.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "ImageIOAnalyze.h"
#include <chrono>

#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

namespace Gadgetron
{

// ------------------------------------------------------------------------
template <typename T>
void save_complex_array_to_file_4d(hoNDArray<std::complex<T>> &array, int dim_1, int dim_2, int dim_3, int dim_4, std::string& filepath) {
    std::string output = "[";
    for (int lll = 0; lll < dim_1; lll++) {
        if (lll > 0) {
            output += ", ";
        }
        output += "[";
        for (int iii = 0; iii < dim_2; iii++) {
            if (iii > 0) {
                output += ", ";
            }
            output += "[";
            for (int jjj = 0; jjj < dim_3; jjj++) {
                if (jjj > 0) {
                    output += ", ";
                }
                output += "[";
                for (int kkk = 0; kkk < dim_4; kkk++) {
                    if (kkk > 0) {
                        output += ", ";
                    }
                    std::complex<T> val = array(lll, iii, jjj, kkk);
                    T rrr = val.real();
                    T mmm = val.imag();
                    std::stringstream stream_rrr;
                    stream_rrr << std::fixed << std::setprecision(50) << rrr;
                    std::string rrr2 = stream_rrr.str();

                    std::stringstream stream_mmm;
                    stream_mmm << std::fixed << std::setprecision(50) << rrr;
                    std::string mmm2 = stream_mmm.str();

                    output += "[";
                    output += rrr2;
                    output += ",";
                    output += mmm2;
                    output += "]";
                }
                output += "]";
            }
            output += "]";
        }
        output += "]";
    }
    output += "]";

    std::ofstream out(filepath);
    out << output;
    out.close();
}

template <typename T>
void save_complex_array_to_file_3d(hoNDArray<std::complex<T>> &array, int dim_1, int dim_2, int dim_3, std::string& filepath) {
    std::string output = "[";
    for (int iii = 0; iii < dim_1; iii++) {
        if (iii > 0) {
            output += ", ";
        }
        output += "[";
        for (int jjj = 0; jjj < dim_2; jjj++) {
            if (jjj > 0) {
                output += ", ";
            }
            output += "[";
            for (int kkk = 0; kkk < dim_3; kkk++) {
                if (kkk > 0) {
                    output += ", ";
                }
                std::complex<T> val = array(iii, jjj, kkk);
                float rrr = val.real();
                float mmm = val.imag();
                std::stringstream stream_rrr;
                stream_rrr << std::fixed << std::setprecision(50) << rrr;
                std::string rrr2 = stream_rrr.str();

                std::stringstream stream_mmm;
                stream_mmm << std::fixed << std::setprecision(50) << mmm;
                std::string mmm2 = stream_mmm.str();

                output += "[";
                output += rrr2;
                output += ",";
                output += mmm2;
                output += "]";
            }
            output += "]";
        }
        output += "]";
    }
    output += "]";

    std::ofstream out(filepath);
    out << output;
    out.close();
}

template EXPORTMRICORE void save_complex_array_to_file_4d(hoNDArray<std::complex<float>> &array, int dim_1, int dim_2, int dim_3, int dim_4, std::basic_string<char>& filepath);
template EXPORTMRICORE void save_complex_array_to_file_3d(hoNDArray<std::complex<float>> &array, int dim_1, int dim_2, int dim_3, std::basic_string<char>& filepath);
template EXPORTMRICORE void save_complex_array_to_file_4d(hoNDArray<std::complex<double>> &array, int dim_1, int dim_2, int dim_3, int dim_4, std::basic_string<char>& filepath);
template EXPORTMRICORE void save_complex_array_to_file_3d(hoNDArray<std::complex<double>> &array, int dim_1, int dim_2, int dim_3, std::basic_string<char>& filepath);


template <typename T> 
void custom_apply_unmix_coeff_aliased_image_3D(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        size_t RO = aliasedIm.get_size(0);
        size_t E1 = aliasedIm.get_size(1);
        size_t E2 = aliasedIm.get_size(2);
        size_t srcCHA = aliasedIm.get_size(3);

        size_t N = aliasedIm.get_size(4);

        std::vector<size_t> dim;
        aliasedIm.get_dimensions(dim);

        GADGET_CHECK_THROW(unmixCoeff.get_size(0) == RO);
        GADGET_CHECK_THROW(unmixCoeff.get_size(1) == E1);
        GADGET_CHECK_THROW(unmixCoeff.get_size(2) == E2);
        GADGET_CHECK_THROW(unmixCoeff.get_size(3) == srcCHA);

        if (complexIm.get_size(0) != RO
            || complexIm.get_size(1) != E1
            || complexIm.get_size(2) != E2
            || complexIm.get_number_of_elements() != RO*E1*E2*N)
        {
            complexIm.create(RO, E1, E2, N);
        }

        hoNDArray<T> buffer;
        buffer.create(dim);

        Gadgetron::multiply(aliasedIm, unmixCoeff, buffer);

        hoNDArray<T> bufferIm;
        bufferIm.create(RO, E1, E2, 1, N, complexIm.begin());


        int RO2 = buffer.get_size(0);
        int E12 = buffer.get_size(1);
        int E22 = buffer.get_size(2);
        int CHA2 = buffer.get_size(3);
        int N2 = buffer.get_size(4);
        int S2 = buffer.get_size(5);
        int SLC2 = buffer.get_size(6);
        GDEBUG_STREAM(">>>>>>>>>>>>>>>>>>>>>>>>>>>> R0: " << RO2 << " | E1: " << E12 << " | E2: " << E22 << " | CHA: " << CHA2 << " | N: " << N2 << " | S: " << S2 << " | SLC: " << SLC2);

        const auto p1 = std::chrono::system_clock::now();
        std::string filepath = "/mnt/i/Data PhD/OneDrive - Imperial College London/PhD/In-house recon/saved_data/temp/ksp_post_grappa_nosms_";
        std::string ext = ".txt";
        std::string time = std::to_string(std::chrono::duration_cast<std::chrono::microseconds>(p1.time_since_epoch()).count());
        std::string fullpath = filepath + time + ext;

        save_complex_array_to_file_4d(buffer, RO2, E12, E22, CHA2, fullpath);

        Gadgetron::sum_over_dimension(buffer, bufferIm, 3);
    }
    catch (...)
    {
        GADGET_THROW("Errors in custom_apply_unmix_coeff_aliased_image_3D(...) ... ");
    }
}

template EXPORTMRICORE void custom_apply_unmix_coeff_aliased_image_3D(const hoNDArray< std::complex<float> >& aliasedIm, const hoNDArray< std::complex<float> >& unmixCoeff, hoNDArray< std::complex<float> >& complexIm);
template EXPORTMRICORE void custom_apply_unmix_coeff_aliased_image_3D(const hoNDArray< std::complex<double> >& aliasedIm, const hoNDArray< std::complex<double> >& unmixCoeff, hoNDArray< std::complex<double> >& complexIm);
}
