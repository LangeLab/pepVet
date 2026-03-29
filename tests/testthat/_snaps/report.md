# BSA trypsin report matches snapshot (Good/Moderate verdict)

    Code
      digest_report(ev)
    Output
      -- sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4  —  tryps
      sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4  v Good  (composite: 0.853)
        S_length   ███████░░░ 0.688
        S_coverage ██████████ 0.997
        S_count    ██████████ 1.000
        S_hydro    ████████░░ 0.769
        S_charge   ████████░░ 0.778
      --------------------------------------------------------------------------------

# Histone H3 trypsin report matches snapshot (Poor verdict)

    Code
      digest_report(ev)
    Output
      -- sp|P68431|H31_HUMAN Histone H3.1 OS=Homo sapiens OX=9606 GN=H3C1 PE=1 SV=2  —
      sp|P68431|H31_HUMAN Histone H3.1 OS=Homo sapiens OX=9606 GN=H3C1 PE=1 SV=2  x Poor  (composite: 0.297)
        S_length   █░░░░░░░░░ 0.133
        S_coverage ██░░░░░░░░ 0.235
        S_count    █░░░░░░░░░ 0.088
        S_hydro    ██████████ 1.000
        S_charge   ██░░░░░░░░ 0.250
      --------------------------------------------------------------------------------

# BSA lysc report matches snapshot

    Code
      digest_report(ev)
    Output
      -- sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4  —  lysc 
      sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4  ! Moderate  (composite: 0.640)
        S_length   ██████░░░░ 0.623
        S_coverage ████████░░ 0.761
        S_count    █████░░░░░ 0.501
        S_hydro    ███████░░░ 0.737
        S_charge   ██████░░░░ 0.553
      --------------------------------------------------------------------------------

# multi-enzyme comparison report matches snapshot

    Code
      digest_report(comp)
    Output
      -- sp|P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4 ---------
        enzyme                            S_len  S_cov  S_cnt  S_hyd  S_chg  composite  verdict
      -----------------------------------------------------------------------------------------
      > glutamyl endopeptidase            0.633  0.827  0.563  0.763  0.947  0.734  Good    
        lysc                              0.623  0.761  0.501  0.737  0.553  0.640  Moderate
        trypsin                           0.532  0.786  0.484  0.762  0.429  0.605  Moderate
        asp-n endopeptidase               0.488  0.445  0.329  0.750  0.950  0.554  Moderate
      --------------------------------------------------------------------------------

# Histone H3 multi-enzyme comparison report matches snapshot

    Code
      digest_report(comp)
    Output
      -- sp|P68431|H31_HUMAN Histone H3.1 OS=Homo sapiens OX=9606 GN=H3C1 PE=1 SV=2 --
        enzyme                            S_len  S_cov  S_cnt  S_hyd  S_chg  composite  verdict
      -----------------------------------------------------------------------------------------
      > asp-n endopeptidase               0.600  0.404  0.375  1.000  1.000  0.626  Moderate
        lysc                              0.429  0.522  0.265  0.833  0.833  0.541  Moderate
        trypsin                           0.133  0.235  0.088  1.000  0.250  0.297  Poor    
      --------------------------------------------------------------------------------

