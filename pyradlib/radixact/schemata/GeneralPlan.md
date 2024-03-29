* PLAN
  * LAST_OPEN_PAGE
  * GENERAL_PLAN
    * PLAN_PROFILE
      * PROCEDURE_TYPE
      * PLAN_NAME
      * CT_TYPE
      * PLAN_STATE
      * PLANNING_ALGORITHM
      * OPERATOR
      * AE_TITLE
      * TREATMENT_MACHINE_TYPE
      * TIMESTAMP
        * DATETIME
    * SYSTEM_INFO
      * SOFTWARE_VERSION
    * PLAN_TEMPLATE
      * ISOCURVESET
        * ISOCURVE
          * DOSE
          * LEVEL
          * COLOR
            * R
            * G
            * B
            * A
          * LINE_WIDTH
          * DISPLAY
    * PATIENT_PROFILE
      * MEDICAL_ID
      * LAST_NAME
      * FIRST_NAME
      * MIDDLE_INITIAL
      * CT_COORD_ORIGIN_X
      * CT_COORD_ORIGIN_Y
      * CT_COORD_ORIGIN_Z
      * CT_COORD_SPACING_X
      * CT_COORD_SPACING_Y
      * CT_COORD_SPACING_Z
      * CT_COORD_SIZE_X
      * CT_COORD_SIZE_Y
      * CT_COORD_SIZE_Z
    * MM_TO_DICOM_TRANSFORM
      * T00
      * T01
      * T02
      * T03
      * T10
      * T11
      * T12
      * T13
      * T20
      * T21
      * T22
      * T23
      * T30
      * T31
      * T32
      * T33
    * IMAGESET
      * IMAGE
        * MODALITY
        * UNIT_TYPE
        * DICOM_IMAGE_STUDY_UID
        * DICOM_IMAGE_SERIES_UID
        * FILEPATH
        * PATIENTID
        * PATIENTLASTNAME
        * WLFILTER_2D
          * WINDOW
        * WLFILTER_3D
        * WLFILTER_DRR
      * DISPLAY_CT_IMAGE_TYPE
    * FIDUCIALSET
      * FIDUCIAL
        * X
        * Y
        * Z
        * FIDUCIAL_AUTOMATICALLY_FOUND
    * VOI_INTERSECTION_SETTINGS
      * VOI
        * ID
        * BEAM_INTERSECTION_CHECK_TYPE
    * AUTOMATED_OPERATIONS
      * IS_QUICK_PLAN
      * FUSION_AUTOMATED
      * DEFINED_ALIGN_CENTER
      * PERFORMED_OPTIMIZATION
      * CALCED_EVAL_DOSE
    * AUTOSEG
      * VOISET
        * VOI
          * MAPPING_TO_AUTO_SEG_STRUCTURE
          * GENERATED_WITH_AUTO_SEGMENTATION
          * CONTOURING_METHOD
    * DOSE_POINTS
      * DOSE_POINTS_SET
    * PLAN_SETUP
      * NUMBER_OF_FRACTIONS
      * PRESCRIPTION_MODE
      * PRESCRIBED_DOSE
      * PRESCRIBED_PERCENTAGE
      * PRESCRIBED_VOI
    * DENSITY_MODEL
      * NAME
      * WARNING1
      * UUID
      * TYPE
      * VERSION
      * EXTENSION
      * DENSITY_TABLE
        * DENSITY_POINT
          * CT
          * ELECTRON_DENSITY
          * MASS_DENSITY
      * WLFILTER_MASS_DENSITY
        * WINDOW
        * LEVEL
      * WLFILTER_ELECTRON_DENSITY
    * DX_VX_VALUES
      * DX_VX_DATA_SET
        * DX_VX_DATA
          * ACTIVE_PLAN_VOI_INDEX
          * SPECIFIED_QUANTITY
          * SPECIFIED_VALUE
          * DX_VX_CRITERIA
            * CRITERIA_QUANTITY
            * CRITERIA_OP
            * CRITERIA_OPERAND1
            * CRITERIA_OPERAND2
            * CRITERIA_WARN_MIN
            * CRITERIA_WARN_MAX