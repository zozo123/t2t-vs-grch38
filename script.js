(function () {
  "use strict";

  // Colour groups in the figures: acrocentric chromosomes and the three
  // chromosomes whose added sequence is mostly satellite heterochromatin.
  var acro = new Set(["chr13", "chr14", "chr15", "chr21", "chr22"]);
  var het = new Set(["chrY", "chr9", "chr16"]);

  function fmt(n) { return Math.round(n).toLocaleString("en-US"); }
  function mb(bp) { return (bp / 1e6).toFixed(1); }
  function one(x) { return x.toFixed(1); }
  // Explicit sign; U+2212 for negatives.
  function signed(x, f) { return (x < 0 ? "−" : "+") + f(Math.abs(x)); }
  function get(obj, path) {
    return path.split(".").reduce(function (o, k) { return o == null ? o : o[k]; }, obj);
  }
  function range(vals, pad) {
    return { lo: Math.min(0, Math.min.apply(null, vals)), hi: Math.max(0, Math.max.apply(null, vals)) * pad };
  }

  function bind(data) {
    var aliases = {
      "totals.net_new_mb": mb(data.totals.net_new_bp) + " Mb",
      "totals.grch38_n_mb": mb(data.totals.grch38_n_bp) + " Mb",
      "totals.extra_islands": fmt(data.totals.extra_islands),
      "totals.new_seq_isl_per_mb": one(data.totals.new_seq_isl_per_mb),
      "timings.warmup_sec": one(data.timings.warmup_sec) + " s",
      "timings.map_wallclock_sec": one(data.timings.map_wallclock_sec) + " s",
      "snapshot.save_sec": one(data.snapshot.save_sec) + " s",
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
    var rows = data.chromosomes;
    var X = range(rows.map(function (r) { return r.delta_bp / 1e6; }), 1.08);
    var Y = range(rows.map(function (r) { return r.delta_islands; }), 1.12);
    host.innerHTML = '<span class="axis y">extra CpG-rule windows</span><span class="axis x">added usable sequence (Mb, net)</span>';
    rows.forEach(function (r) {
      var p = document.createElement("div");
      p.className = "pt" + (acro.has(r.chrom) ? " acro" : het.has(r.chrom) ? " het" : "");
      p.style.left = ((r.delta_bp / 1e6 - X.lo) / (X.hi - X.lo) * 92 + 4) + "%";
      p.style.bottom = ((r.delta_islands - Y.lo) / (Y.hi - Y.lo) * 86 + 6) + "%";
      p.title = r.chrom + ": " + signed(r.delta_bp, mb) + " Mb, " + signed(r.delta_islands, fmt) + " CpG-rule windows";
      p.textContent = r.chrom.replace("chr", "");
      host.appendChild(p);
    });
  }

  // One bar on a track whose zero sits at the left unless the column has losses;
  // losses are drawn leftward from zero with the .neg style.
  function track(v, R, cls, title) {
    var span = R.hi - R.lo, zero = -R.lo / span * 100, w = Math.abs(v) / span * 100;
    var x = v < 0 ? zero - w : zero;
    return '<span class="track" title="' + title + '">' +
      (zero > 0 ? '<span class="zero" style="left:' + zero.toFixed(2) + '%"></span>' : "") +
      '<span class="fill ' + cls + (v < 0 ? " neg" : "") + '" style="--x:' + x.toFixed(2) + "%;--w:" + w.toFixed(2) + '%"></span></span>';
  }
  function val(x, f, unit) {
    var s = signed(x, f) + unit;
    return x < 0 ? '<b class="neg">' + s + "</b>" : s;
  }

  function bars(data) {
    var host = document.getElementById("bars");
    if (!host) return;
    var rows = data.chromosomes.slice().sort(function (a, b) { return b.delta_bp - a.delta_bp; });
    var B = range(rows.map(function (r) { return r.delta_bp; }), 1);
    var I = range(rows.map(function (r) { return r.delta_islands; }), 1);
    host.innerHTML = '<div class="barrow head"><span></span><span>Δ usable sequence</span><span>Δ CpG-rule windows</span><span class="val">Δ Mb · Δ windows</span></div>';
    rows.forEach(function (r) {
      var el = document.createElement("div");
      el.className = "barrow";
      el.innerHTML =
        '<span class="name">' + r.chrom + "</span>" +
        track(r.delta_bp, B, het.has(r.chrom) ? "green" : "", signed(r.delta_bp, mb) + " Mb usable sequence") +
        track(r.delta_islands, I, "hot", signed(r.delta_islands, fmt) + " CpG-rule windows") +
        '<span class="val">' + val(r.delta_bp, mb, " Mb") + " · " + val(r.delta_islands, fmt, "") + "</span>";
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
      el.title = r.chrom + ": " + one(r.grch38.density) + " → " + one(r.t2t.density) + " windows/Mb";
      el.innerHTML =
        '<span class="name">' + r.chrom + "</span>" +
        '<span class="slope" style="--a:' + a + ";--b:" + b + '"><i class="' + (d < 0 ? "down" : "") + '"></i></span>' +
        '<span class="val">' + val(d, one, "/Mb") + "</span>";
      host.appendChild(el);
    });
  }

  fetch("./data/compare.json")
    .then(function (r) {
      if (!r.ok) throw new Error("data/compare.json: HTTP " + r.status);
      return r.json();
    })
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
