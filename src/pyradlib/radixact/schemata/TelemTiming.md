# TelemTiming*.dat

TelemTiming*.dat files are human-readable text files, describing the timing of model build and delivery fragments. The file name typically takes the form:
```
TelemTiming_[SeriesDescription].dat
```
where the SeriesDescription, defined in the DICOM RTRECORD with the SeriesDescription attribute (0x0008,0x103E), has the form:
```
3480-[UID1]-[UID2]-[Fraction]-[Session]
```

The contents of the TelemTiming*.dat files are as shown below, with each line formatted as UTF-8 characters, separated by line breaks:
1. String. Version number. Value is "v3.0" for all observed cases.
2. Integer. Unknown. Value is "1" for all observed cases.
3. Integer. Unknown. Value is "8" for all observed cases.
4. Comma-separated integers, length is number of fragments. Date-time stamps. Value is "0" for all model build fragments. Value is unix time-stamp indicating start of delivery (i.e., seconds since Jan 01 1970 in UTC) for all delivery fragments.
5. Comma-seperated integers, length is number of fragments. Unknown. Value is "0" for all model build fragments. Value is non-"0" for all delivery fragments.
6. Comma-separated integers, length is number of fragments. Start tau or projection. Value is "0" for all model build fragments. Value is starting tau for all delivery fragments.
7. Comma-separated floats, length is number of fragments. End tau or projection. Value is "0" for all model build fragments. Value is ending tau for all delivery fragments.
8. Comma-separated floats, length is number of fragments. Time per projection. Value is "0" for all model build fragments. Value is time, in seconds, for each unit projection for all delivery fragments.
9. Integer. Unknown. Value is starting tau/projection of treatment for all observed cases. Hypothesis is remains constant for interrupted fractions delivered over multiple sessions.
10. Comma-separated floats, length is number of fragments. Values are "0" for all observed cases.
11. Boolean. Unknown. Value is "False" for all observed cases.
