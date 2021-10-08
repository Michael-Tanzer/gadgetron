import numpy as np
import gadgetron
import ismrmrd
import logging
import time
from ismrmrdtools import grappa, transform, coils
from gadgetron.util.cfft import cfftn, cifftn
from gadgetron.types.recon_data import *
# import matplotlib.pyplot as plt
from PIL import Image


def InterceptorGadget1(connection):
    counter = 0
    for acq in connection:
        print(f"----------------------------5  object of class {type(acq)} with dim {acq.data.shape}  5----------------------------")
        # print(f"----------------------------5.1  {acq.meta}  5.1----------------------------")
        # print(f"----------------------------5.2  {acq.headers}  5.2----------------------------")
        for i in range(acq.data.shape[-1]):
            data = acq.data[:, :, 0, 0, 0, 0, i].astype(float)
            print(f"----------------------------6  {data.shape}  6----------------------------")
            im = Image.fromarray(data)
            im = im.convert("L")
            im.save(f"/mnt/i/Data PhD/OneDrive - Imperial College London/PhD/In-house recon/test_save/img_test2_{counter:03d}.png")
            counter += 1
        connection.send(acq)


def InterceptorGadget2(connection: gadgetron.external.connection.Connection):
    hdr = connection.header

    print(f"----------------------------0  {connection.header.toxml()}  0----------------------------")
    print(f"----------------------------0  {[i for i in connection.config.keys()]}  0----------------------------")
    print(f"----------------------------0  {type(connection)}  0----------------------------")

    accf = hdr.encoding[0].parallelImaging.accelerationFactor.kspace_encoding_step_1
    counter = 0
    for acq in connection:
        acq: gadgetron.types.ReconData
        for bit in acq.bits:
            bit: gadgetron.types.recon_data.ReconBit
            print(f"----------------------------1  {bit.data.data.shape}  1----------------------------")
            print(f"----------------------------2  {bit.data.data.dtype}  2----------------------------")
            data = bit.data.data
            print(f"----------------------------3  {data.squeeze().shape}  3----------------------------")
            if len(data.squeeze().shape) == 4: # TODO: this is probs wrong, I need to use the two last dimensions separately for 2 recons and send them
                data = data.squeeze().mean(3).transpose((2, 0, 1))
            else:
                data = data.squeeze().transpose((2, 0, 1))

            # (csm, rho) = coils.calculate_csm_inati_iter(data)
            # (csm, rho) = coils.calculate_csm_walsh(data)
            (unmix, gmap) = grappa.calculate_grappa_unmixing(data, 2, regularization_factor=0.0005)

            recon_data = transform.transform_kspace_to_image(np.squeeze(data), (1, 2)) * np.sqrt(accf)
            recon = np.sum(unmix * recon_data, 0)

            # acq_headers = [item.item() for item in np.nditer(bit.data.headers, ('refs_ok', 'zerosize_ok'), order='F')]
            # meta = connection.meta

            # gadgetron.types.ImageArray(recon, headers=, meta=, waveform=None, acq_headers=acq_headers)

            # im = Image.fromarray(recon.astype(float))
            # im = im.convert("L")
            #
            # img_head = ismrmrd.ImageHeader()
            # img_head.channels = acq.active_channels
            # img_head.slice = acq.idx.slice
            # img_head.matrix_size[0] = hdr.encoding[0].reconSpace.matrixSize.x
            # img_head.matrix_size[1] = hdr.encoding[0].reconSpace.matrixSize.y
            # img_head.matrix_size[2] = hdr.encoding[0].reconSpace.matrixSize.z
            # img_head.position = acq.position
            # img_head.read_dir = acq.read_dir
            # img_head.phase_dir = acq.phase_dir
            # img_head.slice_dir = acq.slice_dir
            # img_head.patient_table_position = acq.patient_table_position
            # img_head.acquisition_time_stamp = acq.acquisition_time_stamp
            # img_head.image_index = self.myCounter
            # img_head.image_series_index = self.mySeries
            # img_head.data_type = ismrmrd.DATATYPE_CXFLOAT

            # im.save(f"/mnt/i/Data PhD/OneDrive - Imperial College London/PhD/In-house recon/test_save/img_test_{counter:04d}.png")
            # counter += 1
            print(f"----------------------------3.5  {recon.shape}  3.5----------------------------")
            connection.send(ismrmrd.Image.from_array(recon))
        # connection.send(acq)


def FilterWrongDimensionGadget(connection):
    for acq in connection:
        print(f"object of class {type(acq)} with dim {acq.data.shape} {'WENT THROUGH' if acq.data.shape[1] == 320 else 'WAS STOPPED'}\n")
        if acq.data.shape[1] == 320:
            connection.send(acq)


def FilterWrongDimensionGadget2(connection):
    hdr = connection.header
    recon_x = hdr.encoding[0].reconSpace.matrixSize.x
    for acq in connection:
        condition = (acq.data.shape[1] == 2*recon_x and isinstance(acq, ismrmrd.acquisition.Acquisition)) or isinstance(acq, ismrmrd.waveform.Waveform)
        if condition:
            connection.send(acq)
        else:
            print(f"object of class {type(acq)} with dim {acq.data.shape} WAS STOPPED")


def Remove2XOversampling(connection):
    for acq in connection:
        orig_size = list(acq.data.shape)
        data2 = acq.data.reshape([acq.data.shape[0], int(acq.data.size / acq.data.shape[0])])
        new_length = data2.shape[1] // 2
        data2 = cfftn(cifftn(data2, axes=[0])[:, (0 + (new_length // 2)): (new_length + (new_length // 2))], axes=[0])
        orig_size[1] = new_length
        data2.reshape(tuple(orig_size))
        acq.resize(*orig_size[::-1])
        acq.center_sample = orig_size[1]
        acq.data[:] = data2
        acq.samples = new_length

        connection.send(acq)

# <?xml version="1.0" ?><ismrmrdHeader xmlns="http://www.ismrm.org/ISMRMRD"><studyInformation><studyTime>17:18:09</studyTime></studyInformation><measurementInformation><measurementID>175749__2.2.1107.5.2.50.175749.30000020113008352930900002402_468</measurementID><patientPosition>HFS</patientPosition><protocolName>SE_b0150500_2SL_SMS2_PAT2_DF100_T6</protocolName><measurementDependency><dependencyType>SenMap</dependencyType><measurementID>175749__2.2.1107.5.2.50.175749.30000020113008352930900002402_432</measurementID></measurementDependency><measurementDependency><dependencyType>Noise</dependencyType><measurementID>175749__2.2.1107.5.2.50.175749.30000020113008352930900002402_432</measurementID></measurementDependency><frameOfReferenceUID>1.3.12.2.1107.5.2.50.175749.1.20201130162006131.0.0.0</frameOfReferenceUID></measurementInformation><acquisitionSystemInformation><systemVendor>Siemens</systemVendor><systemModel>MAGNETOM Vida</systemModel><systemFieldStrength_T>2.89362</systemFieldStrength_T><relativeReceiverNoiseBandwidth>0.793</relativeReceiverNoiseBandwidth><receiverChannels>20</receiverChannels><coilLabel><coilNumber>20</coilNumber><coilName>HeadNeck_20_TCS:1:H12</coilName></coilLabel><coilLabel><coilNumber>19</coilNumber><coilName>HeadNeck_20_TCS:1:H11</coilName></coilLabel><coilLabel><coilNumber>10</coilNumber><coilName>HeadNeck_20_TCS:1:H14</coilName></coilLabel><coilLabel><coilNumber>9</coilNumber><coilName>HeadNeck_20_TCS:1:H13</coilName></coilLabel><coilLabel><coilNumber>2</coilNumber><coilName>HeadNeck_20_TCS:1:H32</coilName></coilLabel><coilLabel><coilNumber>1</coilNumber><coilName>HeadNeck_20_TCS:1:H31</coilName></coilLabel><coilLabel><coilNumber>4</coilNumber><coilName>HeadNeck_20_TCS:1:H34</coilName></coilLabel><coilLabel><coilNumber>3</coilNumber><coilName>HeadNeck_20_TCS:1:H33</coilName></coilLabel><coilLabel><coilNumber>16</coilNumber><coilName>HeadNeck_20_TCS:1:N12</coilName></coilLabel><coilLabel><coilNumber>15</coilNumber><coilName>HeadNeck_20_TCS:1:N11</coilName></coilLabel><coilLabel><coilNumber>18</coilNumber><coilName>HeadNeck_20_TCS:1:H21</coilName></coilLabel><coilLabel><coilNumber>17</coilNumber><coilName>HeadNeck_20_TCS:1:H22</coilName></coilLabel><coilLabel><coilNumber>11</coilNumber><coilName>HeadNeck_20_TCS:1:H23</coilName></coilLabel><coilLabel><coilNumber>12</coilNumber><coilName>HeadNeck_20_TCS:1:H24</coilName></coilLabel><coilLabel><coilNumber>13</coilNumber><coilName>HeadNeck_20_TCS:1:H42</coilName></coilLabel><coilLabel><coilNumber>14</coilNumber><coilName>HeadNeck_20_TCS:1:H41</coilName></coilLabel><coilLabel><coilNumber>6</coilNumber><coilName>HeadNeck_20_TCS:1:H44</coilName></coilLabel><coilLabel><coilNumber>5</coilNumber><coilName>HeadNeck_20_TCS:1:H43</coilName></coilLabel><coilLabel><coilNumber>8</coilNumber><coilName>HeadNeck_20_TCS:1:N22</coilName></coilLabel><coilLabel><coilNumber>7</coilNumber><coilName>HeadNeck_20_TCS:1:N21</coilName></coilLabel><institutionName>Royal Brompton Hospital</institutionName></acquisitionSystemInformation><experimentalConditions><H1resonanceFrequency_Hz>123256020</H1resonanceFrequency_Hz></experimentalConditions><encoding><encodedSpace><matrixSize><x>160</x><y>100</y><z>1</z></matrixSize><fieldOfView_mm><x>320.0</x><y>100.0</y><z>6.0</z></fieldOfView_mm></encodedSpace><reconSpace><matrixSize><x>160</x><y>50</y><z>1</z></matrixSize><fieldOfView_mm><x>320.0</x><y>100.0</y><z>6.0</z></fieldOfView_mm></reconSpace><encodingLimits><kspace_encoding_step_1><minimum>0</minimum><maximum>49</maximum><center>25</center></kspace_encoding_step_1><kspace_encoding_step_2><minimum>0</minimum><maximum>0</maximum><center>0</center></kspace_encoding_step_2><average><minimum>0</minimum><maximum>0</maximum><center>0</center></average><slice><minimum>0</minimum><maximum>1</maximum><center>0</center></slice><contrast><minimum>0</minimum><maximum>0</maximum><center>0</center></contrast><phase><minimum>0</minimum><maximum>0</maximum><center>0</center></phase><repetition><minimum>0</minimum><maximum>0</maximum><center>0</center></repetition><set><minimum>0</minimum><maximum>0</maximum><center>0</center></set><segment><minimum>0</minimum><maximum>0</maximum><center>0</center></segment></encodingLimits><trajectory>epi</trajectory><trajectoryDescription><identifier>ConventionalEPI</identifier><userParameterLong><name>etl</name><value>128</value></userParameterLong><userParameterLong><name>numberOfNavigators</name><value>3</value></userParameterLong><userParameterLong><name>rampUpTime</name><value>130</value></userParameterLong><userParameterLong><name>rampDownTime</name><value>130</value></userParameterLong><userParameterLong><name>flatTopTime</name><value>450</value></userParameterLong><userParameterLong><name>acqDelayTime</name><value>83</value></userParameterLong><userParameterLong><name>numSamples</name><value>320</value></userParameterLong><userParameterDouble><name>dwellTime</name><value>1.7</value></userParameterDouble><comment>Conventional EPI sequence</comment></trajectoryDescription><parallelImaging><accelerationFactor><kspace_encoding_step_1>2</kspace_encoding_step_1><kspace_encoding_step_2>1</kspace_encoding_step_2></accelerationFactor><calibrationMode>other</calibrationMode></parallelImaging></encoding><sequenceParameters><TR>1500.0</TR><TE>45.0</TE><TI>2500.0</TI><sequence_type>EPI</sequence_type><echo_spacing>0.71</echo_spacing></sequenceParameters></ismrmrdHeader>
# <?xml version="1.0"?>\n<ismrmrdMeta>\n\t<meta>\n\t\t<name>GADGETRON_DataRole</name>\n\t\t<value>Image</value>\n\t</meta>\n\t<meta>\n\t\t<name>GADGETRON_ImageComment</name>\n\t\t<value>GT</value>\n\t</meta>\n\t<meta>\n\t\t<name>GADGETRON_ImageNumber</name>\n\t\t<value>1</value>\n\t</meta>\n\t<meta>\n\t\t<name>GADGETRON_ImageProcessingHistory</name>\n\t\t<value>GT</value>\n\t</meta>\n\t<meta>\n\t\t<name>GADGETRON_SeqDescription</name>\n\t\t<value>_GT</value>\n\t</meta>\n\t<meta>\n\t\t<name>PatientPosition</name>\n\t\t<value>3.11422</value>\n\t\t<value>49.6208</value>\n\t\t<value>19.9085</value>\n\t</meta>\n\t<meta>\n\t\t<name>acquisition_time_stamp</name>\n\t\t<value>24915887</value>\n\t</meta>\n\t<meta>\n\t\t<name>encoded_matrix</name>\n\t\t<value>160</value>\n\t\t<value>100</value>\n\t\t<value>1</value>\n\t</meta>\n\t<meta>\n\t\t<name>encoding</name>\n\t\t<value>0</value>\n\t</meta>\n\t<meta>\n\t\t<name>encoding_FOV</name>\n\t\t<value>320</value>\n\t\t<value>100</value>\n\t\t<value>6</value>\n\t</meta>\n\t<meta>\n\t\t<name>gadgetron_sha1</name>\n\t\t<value>1b37cbb8963c7d776d9f3694b312a8c08fda5830</value>\n\t</meta>\n\t<meta>\n\t\t<name>patient_table_position</name>\n\t\t<value>0</value>\n\t\t<value>0</value>\n\t\t<value>-1.119e+06</value>\n\t</meta>\n\t<meta>\n\t\t<name>phase_dir</name>\n\t\t<value>0.318256</value>\n\t\t<value>0.943436</value>\n\t\t<value>0.0929609</value>\n\t</meta>\n\t<meta>\n\t\t<name>physiology_time_stamp</name>\n\t\t<value>1914060</value>\n\t\t<value>0</value>\n\t\t<value>0</value>\n\t</meta>\n\t<meta>\n\t\t<name>read_dir</name>\n\t\t<value>-0.94732</value>\n\t\t<value>0.32022</value>\n\t\t<value>-0.00663849</value>\n\t</meta>\n\t<meta>\n\t\t<name>recon_FOV</name>\n\t\t<value>320</value>\n\t\t<value>100</value>\n\t\t<value>6</value>\n\t</meta>\n\t<meta>\n\t\t<name>recon_matrix</name>\n\t\t<value>160</value>\n\t\t<value>50</value>\n\t\t<value>1</value>\n\t</meta>\n\t<meta>\n\t\t<name>sampling_limits_E1</name>\n\t\t<value>25</value>\n\t\t<value>50</value>\n\t\t<value>74</value>\n\t</meta>\n\t<meta>\n\t\t<name>sampling_limits_E2</name>\n\t\t<value>0</value>\n\t\t<value>0</value>\n\t\t<value>0</value>\n\t</meta>\n\t<meta>\n\t\t<name>sampling_limits_RO</name>\n\t\t<value>0</value>\n\t\t<value>80</value>\n\t\t<value>159</value>\n\t</meta>\n\t<meta>\n\t\t<name>slice_dir</name>\n\t\t<value>-0.036031</value>\n\t\t<value>-0.085951</value>\n\t\t<value>0.995648</value>\n\t</meta>\n</ismrmrdMeta>\n', '<?xml version="1.0"?>\n<ismrmrdMeta>\n\t<meta>\n\t\t<name>GADGETRON_DataRole</name>\n\t\t<value>Image</value>\n\t</meta>\n\t<meta>\n\t\t<name>GADGETRON_ImageComment</name>\n\t\t<value>GT</value>\n\t</meta>\n\t<meta>\n\t\t<name>GADGETRON_ImageNumber</name>\n\t\t<value>2</value>\n\t</meta>\n\t<meta>\n\t\t<name>GADGETRON_ImageProcessingHistory</name>\n\t\t<value>GT</value>\n\t</meta>\n\t<meta>\n\t\t<name>GADGETRON_SeqDescription</name>\n\t\t<value>_GT</value>\n\t</meta>\n\t<meta>\n\t\t<name>PatientPosition</name>\n\t\t<value>3.54659</value>\n\t\t<value>50.6522</value>\n\t\t<value>7.96068</value>\n\t</meta>\n\t<meta>\n\t\t<name>acquisition_time_stamp</name>\n\t\t<value>24916679</value>\n\t</meta>\n\t<meta>\n\t\t<name>encoded_matrix</name>\n\t\t<value>160</value>\n\t\t<value>100</value>\n\t\t<value>1</value>\n\t</meta>\n\t<meta>\n\t\t<name>encoding</name>\n\t\t<value>0</value>\n\t</meta>\n\t<meta>\n\t\t<name>encoding_FOV</name>\n\t\t<value>320</value>\n\t\t<value>100</value>\n\t\t<value>6</value>\n\t</meta>\n\t<meta>\n\t\t<name>gadgetron_sha1</name>\n\t\t<value>1b37cbb8963c7d776d9f3694b312a8c08fda5830</value>\n\t</meta>\n\t<meta>\n\t\t<name>patient_table_position</name>\n\t\t<value>0</value>\n\t\t<value>0</value>\n\t\t<value>-1.119e+06</value>\n\t</meta>\n\t<meta>\n\t\t<name>phase_dir</name>\n\t\t<value>0.318256</value>\n\t\t<value>0.943436</value>\n\t\t<value>0.0929609</value>\n\t</meta>\n\t<meta>\n\t\t<name>physiology_time_stamp</name>\n\t\t<value>1914852</value>\n\t\t<value>0</value>\n\t\t<value>0</value>\n\t</meta>\n\t<meta>\n\t\t<name>read_dir</name>\n\t\t<value>-0.94732</value>\n\t\t<value>0.32022</value>\n\t\t<value>-0.00663849</value>\n\t</meta>\n\t<meta>\n\t\t<name>recon_FOV</name>\n\t\t<value>320</value>\n\t\t<value>100</value>\n\t\t<value>6</value>\n\t</meta>\n\t<meta>\n\t\t<name>recon_matrix</name>\n\t\t<value>160</value>\n\t\t<value>50</value>\n\t\t<value>1</value>\n\t</meta>\n\t<meta>\n\t\t<name>sampling_limits_E1</name>\n\t\t<value>25</value>\n\t\t<value>50</value>\n\t\t<value>74</value>\n\t</meta>\n\t<meta>\n\t\t<name>sampling_limits_E2</name>\n\t\t<value>0</value>\n\t\t<value>0</value>\n\t\t<value>0</value>\n\t</meta>\n\t<meta>\n\t\t<name>sampling_limits_RO</name>\n\t\t<value>0</value>\n\t\t<value>80</value>\n\t\t<value>159</value>\n\t</meta>\n\t<meta>\n\t\t<name>slice_dir</name>\n\t\t<value>-0.036031</value>\n\t\t<value>-0.085951</value>\n\t\t<value>0.995648</value>\n\t</meta>\n</ismrmrdMeta>\n


class TestClassGadget(gadgetron.Gadget):
    def __init__(self,next_gadget=None):
        super(TestClassGadget,self).__init__(next_gadget=next_gadget)
        self.counter = []

    def process_config(self, cfg):
        pass

    def process(self, acq, data, *args):
        print("hello from the most amazing gadget ever\n"*100)