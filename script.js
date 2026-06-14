(function () {
  "use strict";

  var acro = new Set(["chr13", "chr14", "chr15", "chr21", "chr22"]);
  var het = new Set(["chrY", "chr9"]);

  function fmt(n) { return Math.round(n).toLocaleString("en-US"); }
  function mb(bp) { return (bp / 1e6).toFixed(1); }
  function get(obj, path) {
    return path.split(".").reduce(function (o, k) { return o == null ? o : o[k]; }, obj);
  }
  function bind(data) {
    var aliases = {
      "totals.net_new_mb": mb(data.totals.net_new_bp) + " Mb",
      "totals.grch38_n_mb": mb(data.totals.grch38_n_bp) + " Mb",
      "totals.extra_islands": fmt(data.totals.extra_islands),
      "totals.new_seq_isl_per_mb": data.totals.new_seq_isl_per_mb + "/Mb",
      "timings.warmup_sec": data.timings.warmup_sec.toFixed(1) + " s",
      "timings.map_wallclock_sec": data.timings.map_wallclock_sec.toFixed(1) + " s",
      "snapshot.save_sec": data.snapshot.save_sec.toFixed(1) + " s",
      "snapshot.size": data.snapshot.size
    };
    document.querySelectorAll("[data-bind]").forEach(function (el) {
      var k = el.getAttribute("data-bind");
      var v = Object.prototype.hasOwnProperty.call(aliases, k) ? aliases[k] : get(data, k);
      if (v != null) el.textContent = v;
    });
  }

  function scatter(data) {
    var host = document.getElementById("scatter");
    if (!host) return;
    var rows = data.chromosomes.slice();
    var maxX = Math.max.apply(null, rows.map(function (r) { return Math.max(0, r.delta_bp / 1e6); })) * 1.08;
    var maxY = Math.max.apply(null, rows.map(function (r) { return Math.max(0, r.delta_islands); })) * 1.12;
    host.innerHTML = '<span class="axis y">extra candidate islands</span><span class="axis x">new usable sequence (Mb)</span>';
    rows.forEach(function (r) {
      var x = Math.max(0, r.delta_bp / 1e6) / maxX * 92 + 4;
      var y = Math.max(0, r.delta_islands) / maxY * 86 + 6;
      var p = document.createElement("div");
      p.className = "pt" + (acro.has(r.chrom) ? " acro" : het.has(r.chrom) ? " het" : "");
      p.style.left = x + "%";
      p.style.bottom = y + "%";
      p.title = r.chrom + ": +" + mb(r.delta_bp) + " Mb, " + fmt(r.delta_islands) + " islands";
      p.textContent = r.chrom.replace("chr", "");
      host.appendChild(p);
    });
  }

  function bars(data) {
    var host = document.getElementById("bars");
    if (!host) return;
    var rows = data.chromosomes.slice().sort(function (a, b) {
      return b.delta_bp - a.delta_bp;
    });
    var maxBp = Math.max.apply(null, rows.map(function (r) { return Math.abs(r.delta_bp); }));
    var maxIsl = Math.max.apply(null, rows.map(function (r) { return Math.abs(r.delta_islands); }));
    host.innerHTML = "";
    rows.forEach(function (r) {
      var el = document.createElement("div");
      el.className = "barrow";
      var bpW = Math.max(0, Math.abs(r.delta_bp) / maxBp * 100).toFixed(1) + "%";
      var islW = Math.max(0, Math.abs(r.delta_islands) / maxIsl * 100).toFixed(1) + "%";
      el.innerHTML =
        '<span class="name">' + r.chrom + '</span>' +
        '<span class="track" title="new sequence"><span class="fill ' + (het.has(r.chrom) ? "green" : "") + '" style="--w:' + bpW + '"></span></span>' +
        '<span class="track" title="extra islands"><span class="fill hot" style="--w:' + islW + '"></span></span>' +
        '<span class="val">+' + mb(r.delta_bp) + ' Mb</span>';
      host.appendChild(el);
    });
  }

  function density(data) {
    var host = document.getElementById("density");
    if (!host) return;
    var rows = data.chromosomes.slice().sort(function (a, b) {
      return (b.t2t.density - b.grch38.density) - (a.t2t.density - a.grch38.density);
    });
    var maxD = Math.max.apply(null, rows.map(function (r) { return Math.max(r.grch38.density, r.t2t.density); }));
    host.innerHTML = "";
    rows.forEach(function (r) {
      var a = (r.grch38.density / maxD * 94 + 2).toFixed(1) + "%";
      var b = (r.t2t.density / maxD * 94 + 2).toFixed(1) + "%";
      var d = r.t2t.density - r.grch38.density;
      var el = document.createElement("div");
      el.className = "density-row";
      el.innerHTML =
        '<span class="name">' + r.chrom + '</span>' +
        '<span class="slope" style="--a:' + a + ';--b:' + b + '"><i class="' + (d < 0 ? "down" : "") + '"></i></span>' +
        '<span class="val">' + (d >= 0 ? "+" : "") + d.toFixed(1) + '/Mb</span>';
      host.appendChild(el);
    });
  }

  fetch("./data/compare.json")
    .then(function (r) { return r.json(); })
    .then(function (data) {
      bind(data);
      scatter(data);
      bars(data);
      density(data);
    })
    .catch(function (err) {
      console.error(err);
      document.body.classList.add("data-error");
    });
})();
