# BSA trypsin report matches snapshot (Good/Moderate verdict)

    Code
      digest_report(ev)
    Output
      pepVet digest check
      -------------------
      Protein            sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB
                         PE=1 SV=4
      Enzyme             trypsin
      Preset             standard
      Missed cleavages   Up to 1
      Peptides           157 total; 108 within 7-25 aa
      Verdict            Good
      Composite          0.885
      Component        Score  Profile
      S_length         0.688  [#######---]
      S_coverage       0.997  [##########]
      S_count          1.000  [##########]
      S_hydro          0.769  [########--]
      S_charge         0.778  [########--]

# Histone H3 trypsin report matches snapshot (Moderate verdict)

    Code
      digest_report(ev)
    Output
      pepVet digest check
      -------------------
      Protein            sp|P68431|H31_HUMAN Histone H3.1 OS=Homo sapiens OX=9606
                         GN=H3C1 PE=1 SV=2
      Enzyme             trypsin
      Preset             standard
      Missed cleavages   Up to 1
      Peptides           59 total; 18 within 7-25 aa
      Verdict            Moderate
      Composite          0.619
      Component        Score  Profile
      S_length         0.305  [###-------]
      S_coverage       0.632  [######----]
      S_count          0.662  [#######---]
      S_hydro          0.833  [########--]
      S_charge         0.833  [########--]

# BSA lysc report matches snapshot

    Code
      digest_report(ev)
    Output
      pepVet digest check
      -------------------
      Protein            sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB
                         PE=1 SV=4
      Enzyme             lysc
      Preset             standard
      Missed cleavages   Up to 1
      Peptides           121 total; 76 within 7-25 aa
      Verdict            Good
      Composite          0.823
      Component        Score  Profile
      S_length         0.628  [######----]
      S_coverage       0.857  [#########-]
      S_count          1.000  [##########]
      S_hydro          0.763  [########--]
      S_charge         0.776  [########--]

# multi-enzyme comparison report matches snapshot

    Code
      digest_report(comp)
    Output
      pepVet enzyme comparison
      ------------------------
      Protein            sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB
                         PE=1 SV=4
      Best score         trypsin (0.885, Good)
      Rank Enzyme                 S_len S_cov S_cnt S_hyd S_chg Score Verdict
      -----------------------------------------------------------------------
         1 trypsin                0.688 0.997 1.000 0.769 0.778 0.885 Good
         2 glutamyl endopeptidase 0.681 0.936 1.000 0.802 0.975 0.884 Good
         3 lysc                   0.628 0.857 1.000 0.763 0.776 0.823 Good
         4 asp-n endopeptidase    0.494 0.529 0.923 0.725 0.925 0.673 Good

# Histone H3 multi-enzyme comparison report matches snapshot

    Code
      digest_report(comp)
    Output
      pepVet enzyme comparison
      ------------------------
      Protein            sp|P68431|H31_HUMAN Histone H3.1 OS=Homo sapiens OX=9606
                         GN=H3C1 PE=1 SV=2
      Best score         lysc (0.769, Good)
      Rank Enzyme              S_len S_cov S_cnt S_hyd S_chg Score Verdict
      --------------------------------------------------------------------
         1 lysc                0.593 0.735 1.000 0.625 0.938 0.769 Good
         2 trypsin             0.305 0.632 0.662 0.833 0.833 0.619 Moderate
         3 asp-n endopeptidase 0.333 0.404 0.640 1.000 1.000 0.578 Moderate

