# *.dplan

FinalDeliveryPlan*.dplan and TelemFluence*.dplan files planned and telemetry beam delivery data. The files are human-readable, in that they are ASCII encoded, however they contain a large block of binary-to-ASCII hexadecimal encoded data. The structure is similar to that of INI or TOML files, containing key=value pairs.

The contents of the TelemTiming*.dat files are as shown below, with lines formatted as UTF-8 characters, separated by line breaks:
1. Version key and value. "version=0" for all observed cases.
2. Expected use key and value. "expectedUse=treatment" for all observed cases. Hypothesise that quality assurance could be another expected use.
3. Couch offset reference position key and value. "couchOffsetReferencePosition=isocenter" for all observed cases.
4. Delimiter indicating start of delivery fragment.
5. Seconds per tau key and value for fragment, where tau is a measure of projections.
6. Start tau key and value for fragment.
7. End tau key and value for fragment.
8. Count key and value for fragment, where the the first 64 integer count values indicate the instances of non-zero leaf-open times for each of the 64 leaves, followed by 7 additional values.
9. Data key and value for fragment. Data is encoded using binary-to-ASCII hexadecimal encoding, delimitered with "\<\<EOB\>\>" and "EOB" strings. Data is contained over multiple lines, each containing 64 hexadecimal characters (32 bytes). The structure of this data is summarised later in this file.
10. Delimiter indicating start of planning metadata.
11. Planning metadata fragment ID key and value.
12. Planning metadata counts key and value. Value is array of 0 values for all observed cases.
13. Planning metadata data key and value. Data is encoded using binary-to-ASCII hexadecimal encoding, delimitered with "\<\<EOB\>\>" and "EOB" strings. No data has been present for all observed cases.
14. Delimiter indicating end of planning metadata. 
15. Delimiter indicating start of delivery metadata.
16. Trigger count key and value for delivery. 
17. Delimiter indicating end of delivery metadata.
18. Delimiter indicating end of delivery fragment.

The encoded fragment data is an array of float32 values, with 7 distinct blocks of data:
1. 64 arrays of paired leaf-open and leaf-close tau values, effectively providing leaf-open times for sinogram production. Unlike CSV or image exports, this data is zero-skipped, that is, leaf positions are not defined for projections where leaves remain closed (and so some arrays have no values). The length of each array is defined by the channel count earlier in the dplan file. The total number of leaf-open and leaf-close tau value pairs is equal to the sum of the channel counts for each of the 64 leaf channels.
2. 1 array containing (float32-encoded) tau values in sequential range from 0 to the first non-channel count value (i.e., the 65th channel, maximum tau). These tau values extend beyond the start and end tau values above.
3. 1 array of tau and IEC-X couch position values. There are two paired values, for tau 0 and maximum tau. The start and end IEC-X couch position value has been the same for all observed cases.
4. 1 array of tau and IEC-Y couch position values. There are two paired values, for tau 0 and maximum tau.
5. 1 array of tau and IEC-Z couch position values. There are two paired values, for tau 0 and maximum tau. The start and end IEC-Z couch position value has been the same for all observed cases. There are two paired values, for tau 0 and maximum tau. 
6. 1 array of tau and gantry angle values. The starting gantry angle is -51/180 or approximately -3.5, i.e., the starting angle for a projection centred on gantry angle 0, for all observed cases. The maximum gantry angle is a cumulative measure of total rotation, i.e., it exceeds 360 degrees.
7. 2 arrays of tau and jaw position values, for both jaws, indicating position of jaws during treatment delivery. These are not correctly defined for Synchrony treatments.