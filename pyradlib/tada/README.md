# Treatment and Dose Assessor
This subpackage provides functionality for the assessment of treatment dosimetric quality and deliverability. The Treatment and Dose Assessor (TADA) code was previously implemented in C# by Scott Crowe during post-doctoral research fellowship at the Queensland University of Technology.

## References
The application of the TADA code has been described in:
> S. B. Crowe, T. Kairn, N. Middlebrook, B. Hill, D. R. H. Christie, R. T. Knight, J. Kenny, C. M. Langton, J. V. Trapp, "Retrospective evaluation of dosimetric quality for prostate carcinomas treated with 3D conformal, intensity modulated and volumetric modulated arc radiotherapy," J. Med. Radiat. Sci. 60(4): 131-138 (2013).

> S. B. Crowe, T. Kairn, J. Kenny, R. T. Knight, B. Hill, C. M. Langton, J. V. Trapp, "Treatmnt plan complexity metrics for predicting IMRT pre-treatment quality assurance results," Australas. Phys. Eng. Sci. Med. 37(3): 475-482 (2014).

> S. B. Crowe, T. Kairn, N. Middlebrook, B. Sutherland, B. Hill, J. Kenny, C. M. Langton, J. V. Trapp, "Examination of the properties of IMRT and VMAT beams and evaluation against pre-treatment quality assurance results," Phys. Med. Biol. 60(6): 2587-2601 (2015).

> T. Kairn, D. Papworth, S. B. Crowe, J. Anderson, D. R. H. Christie, "Dosimetric quality, accuracy, and deliverability of modulated radiotherapy treatments for spinal metastases," Med. Dosim. 41(3): 258-266 (2016).

> T. Kairn, S. B. Crowe, C. M. Langton, J. V. Trapp, "Bulk evaluation and comparison of radiotherapy treatment plans for breast cancer," Australas. Phys. Eng. Sci. Med. 39(3): 633-644 (2016).

### Complexity metrics
Many of the complexity metrics implemented in this code have been described in review articles:
> M. Antoine, F. Ralite, C. Soustiel, T. Marsac, P. Sargos, A. Cugny, J. Caron, "Use of metrics to quantify IMRT and VMAT treatment plan complexity: A systematic review and perspectives," Physica Medical 64: 98-108 (2019)

> S. Chiavassa, I. Bessieres, M. Edouard, M. Mathot, A. Moignier, "Complexity metrics for IMRT and VMAT plans: a review of current literature and applications," Br. J. Radiol. 92(1102): 20190270 (2019)

The average leaf pair opening (ALPO) is described in:
> P. Zygmanski, J. H. Kung, "Method of identifying dynamic multileaf collimator irradiation that is highly sensitive to a systematic MLC calibration error," Med. Phys. 28(11): 2220-2226 (2001)

The modulation complexity score (MCS), leaf sequence variability (LSV) and aperture area variability (AAV) are described in:
> A. L. McNiven, M. B. Sharpe, T. G. Purdie, "A new metric for assessing IMRT modulation complexity and plan deliverability," Med. Phys. 37(2): 505-515 (2010)

The mean asymmetry distance (MAD), mean field area (MFA), closed leaf score (CLS), cross-axis score (CAS) and small aperture score (SAS) are described in:
> S. B. Crowe, T. Kairn, J. Kenny, R. T. Knight, B. Hill, C. M. Langton, J. V. Trapp, "Treatment plan complexity metrics for predicting IMRT pre-treatment quality assurance results," Australas. Phys. Eng. Sci. Med. 37(3): 475-482 (2014)

### Conformity indices
Many of the conformity indices implemented in this code were described in the Feuvret review:
> L. Feuvret, G. Noel, J-J. Mazeron, P. Bey, "Conformity index: a review," Int. J. Radiat. Oncol. Biol. Phys. 64(2): 333-342 (2006).

The conformation number (CN), also known as Paddick's conformity index (PCI), is described in:
> A. van't Riet, A. C. Mak, M. A. Moerland, L. H. Elders, W. A. van der Zee, "A conformation number to quantify the degree of conformality in brachytherapy and external beam irradiation: application to the prostate," Int. J. Radiat. Oncol. Biol. Phys. 73(3): 731-736 (1997)

