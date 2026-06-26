# Radixact
This subpackage provides functionality for the analysis of data produced by the Radixact treatment delivery system, including data archived using the Patient Data Extractor software, or cached within the Delivery Analysis software, both provided by Accuray.

The following terms are used in describing granularity of Radixact delivery data:
- Patient.
- Plan.
- Fraction.
- Session, describing an instance of a fraction being opened for delivery at the treatment console. During a session, the fractional dose may be delivered completely, partially, or not at all.
- Delivery segment, to describe uninterrupted use of the system within a session. Delivery segments are delineated by user pauses or recoverable errors in delivery.
- Period, for Synchrony treatments, to distinguish between time building a correlation model and time delivering the treatment, within a given delivery segment.

## Displacement metrics and figures
The plan complexity metrics have been taken from or otherwise inspired by the articles below.
> Cavinato S, Fusella M, Paiusco M, Scaggion A (2023). Quantitative assessment of helical tomotherapy plans complexity. *J Appl Clin Med Phys*, 24(1), e13781. DOI:[10.1002/acm2.13781](https://doi.org/10.1002/acm2.13781).
> 

The displacement distribution metrics and figures implemented have been taken from or otherwise inspired by the articles below.
R<sub>80</sub>, R<sub>90</sub>, R<sub>95</sub>.
> Li HS, Chetty IJ, Enke CA, Foster RD, Willoughby TR, Kupellian PA, Solberg TD (2008). Dosimetric consequences of intrafraction prostate motion. *Int J Radiat Oncol Bio Phys*, 71(3), 801-812. DOI:[10.1016/j.ijrobp.2007.10.049](https://doi.org/10.1016/j.ijrobp.2007.10.049).