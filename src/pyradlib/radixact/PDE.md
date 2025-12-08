# Patient data extractor archive data
Patient Data Extractor archives are a collection of .cdms files consisting of data relating to the patient plan and treatment. These .cdms files, which have FileToken names, are each associated with an original file located within the Accuray environment. This markdown document attempts to highlight where and what these files are.

## Path tokens
In describing directories and filenames below, a number of tokens are used for variable paths.

``[YYMMDD-hhmmss]`` represents a date-time stamp corresponding to the relevant data, e.g., ``20250101-120000`` for midday on January 1st 2025.

``[SeriesInstanceUID]`` represents a DICOM Series Instance Unique Identifier. The UID will be composed of the following, separated by ``.``:
 - ``1.2.840.114358``, an identifier specific to Accuray.
 - a unique patient identifier
 - a date-time stamp ``YYYYMMDD-hhmmss``, as above.
 - a session indicator, e.g., ``1``

``[Slice]`` represents the slice number of a transverse slice of a CT data set.

## Treatment Delivery Console (TDC) data
Subdirectories ordinarily contained in the ``C:\accuray\persistence\tdc\`` directory include:
 - ``encounters\[SeriesInstanceUID]``, with treatment encounter data.
 - ``session\data\Scan_[SeriesInstanceUID]\ImageSeries\``, with daily on-board imaging data.
 - ``session\data\[SeriesInstanceUID]\``, with Synchrony tracking data.
 - ``session\data\CompressedDetectorData\[SeriesInstanceUID]``.
 - ``session\workflow\data\[YYYYMMDD-hhmmss]`` with supporting machine data for treatment sessions.

Files in the ``session\data\Scan_[SeriesInstanceUID]\ImageSeries\`` subdirectory include:
 - ``[SeriesInstanceUID].[Slice]`` DICOM CT image data.

Files in the ``session\data\[SeriesInstanceUID]\`` subdirectory include:
 - ``MotionData.xml`` describing target displacement and Synchrony metrics. Schemata provided in ``motionData.md`` file.
 
Files in the ``session\data\CompressedDetectorData\[SeriesInstanceUID]`` subdirectory include:
 - ``CompressedDetector_X.dat``.
 - ``DetectorMetaData_X.dat``.

Files in the ``session\workflow\data\[YYYYMMDD-hhmmss]`` subdirectory may include:
 - ``airscan.header`` and ``airscan.img``, containing the air scan.
 - ``beam_checksum``.
 - ``beamTriggerConfig.csv``.
 - ``checksum.settings`` and ``checksum.manifest``.
 - ``cone.header`` and ``cone.img``.
 - ``configuration.config``, containing various machine data including imaging calibration.
 - ``controls.config``.
 - ``cor.config``, containing centre of rotation offset.
 - ``couchImage.header`` and ``couchImage.img``.
 - ``dcsCalibrationParams.xml`` and ``dcsControlLoop.xml``, containing machine output calibration data.
 - ``dcom-all.config``, ``dcom-twinnable.config`` and ``dcom-nontwinnable.config``.
 - ``fat.header`` and ``fat.img``.
 - ``files.manifest``.
 - ``galilConfig.csv``.
 - ``jaw2field.config``.
 - ``jfof.config``.
 - ``leafLatency.config``.
 - ``lft.header`` and ``lft.img``.
 - ``measurements.config``.
 - ``parametricKernel.config``.
 - ``penumbra.header`` and ``penumbra.img``.
 - ``shared.config``.
 - ``spectralData.header`` and ``spectralData.img``.
 - ``tagp.config``.
 - ``tomoCouchSettings.config``.
 - ``tomoMachineSettings.config``.

## Clinical Data Management System (CDMS) data
Subdirectories ordinarily included in the ``C:\cdms\`` directory include:
 - ``dicom\STANDARD_CT\``, containing DICOM CT data used for planning.
 - ``dicom\STANDARD_RT_STRUCTURE_SET\``, containing DICOM RTSTRUCT data used for planning.
 - ``tmp\``, containing DICOM REG image registrations and DICOM RTPLAN planning data.

## Delivery Analysis (DA) data
Files that may be in the ``C:\tomo\da\pts\*\`` directory include:
 - ``DetectorByProj_*.det``, containing binary float32 ``>f4`` data of exit detector sinogram.
 - ``TelemFluence_*.dplan``, containing and binary-to-ASCII float64 ``>f8`` data with machine telemetry data.
 - ``TelemTiming_*.dat``

These files will only be present if the patient treatment data has been opened in the Delivery Analysis software prior to use of the PDE tool.

## Temporary data
Files tht may be in the ``C:\accuray\tmp\1\plan\*\`` directory include:
 - ``AutoSegData~*.xml``.
 - ``CouchReplacedDataset\[SeriesInstanceUID].[Slice]``, containing DICOM CT images used for treatment planning.
 - ``DoseData~*.xml``.
 - ``FinalDeliveryPlan~*.dplan``.
 - ``GeneralPlan~*.xml``.
 - ``OptimizationDeliveryPlan~*.dplan``.
 - ``reports.xml``.
 - ``rtss~*.dcm``.
 - ``Sinogram0.png``.
 - ``StaticCouchDeliveryPlan~*.dplan``.
 - ``TomoFinalDoseInfo_rtdose~*.dcm``.
 - ``TomoPlanDetails~*.xml``.
 - ``VOISet~*.xml``.
 - ``VoiSetWithOrthogonalContours~*.xml``.