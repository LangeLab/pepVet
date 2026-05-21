document.addEventListener("DOMContentLoaded", function () {
  var refs = [
    "annotate_cleavage_sites", "digest_protein",
    "score_peptides", "pepvet_preset",
    "calculate_peptide_mass", "calculate_pI",
    "evaluate_digest", "compare_digests", "recommend_enzyme",
    "batch_evaluate", "batch_compare_enzymes", "summarize_batch", "triage_proteins",
    "plot_digest_profile", "plot_coverage_map", "plot_cleavage_map",
    "plot_peptide_overlap_map", "plot_length_distribution", "plot_gravy_landscape",
    "plot_pI_distribution", "plot_mz_distribution", "plot_missed_cleavage_impact",
    "plot_enzyme_comparison", "plot_proteome_overview", "plot_batch_comparison",
    "pepvet_plot_config", "pepvet_plot_config_reset", "pepvet_theme_manuscript",
    "pepvet_theme_presentation", "pepvet_save_figure",
    "digest_report", "pepvet_check",
    "export_peptide_list",
    "aa_properties",
    "pepVet-package"
  ];

  var path = window.location.pathname;
  if (!path.includes("/reference/")) return;
  var match = path.match(/\/([^/]+)\.html$/);
  if (!match) return;
  var current = match[1];
  var idx = refs.indexOf(current);
  if (idx === -1) return;

  var mainEl = document.querySelector("main#main");
  if (!mainEl) return;

  var prev = idx > 0 ? refs[idx - 1] : null;
  var next = idx < refs.length - 1 ? refs[idx + 1] : null;

  var navDiv = document.createElement("div");
  navDiv.className = "ref-prevnext";
  navDiv.style.cssText = "display:flex;gap:0;margin-top:2rem";

  if (prev) {
    var a = document.createElement("a");
    a.href = prev + ".html";
    a.className = "ref-nav-link";
    a.style.cssText = "display:flex;flex-direction:column;width:50%";
    if (!next) a.style.cssText += ";margin-left:auto";
    var dir = document.createElement("span");
    dir.className = "ref-nav-dir";
    dir.textContent = "← Previous";
    var fn = document.createElement("span");
    fn.className = "ref-nav-fn";
    var code = document.createElement("code");
    code.textContent = prev + "()";
    fn.appendChild(code);
    a.appendChild(dir);
    a.appendChild(fn);
    navDiv.appendChild(a);
  }

  if (next) {
    var a = document.createElement("a");
    a.href = next + ".html";
    a.className = "ref-nav-link";
    a.style.cssText = "display:flex;flex-direction:column;width:50%;margin-left:auto;text-align:right";
    var dir = document.createElement("span");
    dir.className = "ref-nav-dir";
    dir.textContent = "Next →";
    var fn = document.createElement("span");
    fn.className = "ref-nav-fn";
    var code = document.createElement("code");
    code.textContent = next + "()";
    fn.appendChild(code);
    a.appendChild(dir);
    a.appendChild(fn);
    navDiv.appendChild(a);
  }

  mainEl.appendChild(navDiv);
});
