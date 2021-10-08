
/** \file   mri_core_grappa.h
    \brief  GRAPPA implementation for 2D and 3D MRI parallel imaging
    \author Hui Xue
*/

#pragma once

#include "mri_core_export.h"
#include "hoNDArray.h"

namespace Gadgetron {
    /// aliasedIm : [RO E1 E2 srcCHA ...]
    template <typename T> EXPORTMRICORE void save_complex_array_to_file_4d(hoNDArray<std::complex<T>> &array, int dim_1, int dim_2, int dim_3, int dim_4, std::string& filepath);
    template <typename T> EXPORTMRICORE void save_complex_array_to_file_3d(hoNDArray<std::complex<T>> &array, int dim_1, int dim_2, int dim_3, std::string& filepath);
    template <typename T> EXPORTMRICORE void custom_apply_unmix_coeff_aliased_image_3D(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);

}
