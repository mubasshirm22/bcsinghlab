/* lab.js — Singh Lab site interactions */

(function () {
  'use strict';

  // ── Mobile nav toggle ──
  var hamburger = document.getElementById('hamburger');
  var nav = document.getElementById('lab-nav');
  if (hamburger && nav) {
    hamburger.addEventListener('click', function () {
      nav.classList.toggle('open');
    });
    // Close when clicking a link
    nav.querySelectorAll('a').forEach(function (link) {
      link.addEventListener('click', function () {
        nav.classList.remove('open');
      });
    });
  }

  // ── Active nav link ──
  var path = window.location.pathname;
  var navLinks = document.querySelectorAll('.lab-nav a');
  navLinks.forEach(function (link) {
    var href = link.getAttribute('href');
    if (href && href !== '/' && path.startsWith(href)) {
      link.classList.add('active');
    } else if (href === '/' && path === '/') {
      link.classList.add('active');
    }
  });

  // ── Header: frosted condense + scroll progress bar ──
  var header = document.getElementById('lab-header');
  var progressBar = document.getElementById('scroll-progress-bar');
  var ticking = false;

  function updateHeader() {
    var scrollY = window.pageYOffset || document.documentElement.scrollTop || 0;
    if (header) {
      if (scrollY > 12) header.classList.add('scrolled');
      else header.classList.remove('scrolled');
    }
    if (progressBar) {
      var doc = document.documentElement;
      var max = (doc.scrollHeight - doc.clientHeight);
      var pct = max > 0 ? (scrollY / max) * 100 : 0;
      progressBar.style.width = Math.min(100, Math.max(0, pct)) + '%';
    }
    ticking = false;
  }

  function onScroll() {
    if (!ticking) {
      window.requestAnimationFrame(updateHeader);
      ticking = true;
    }
  }

  if (header || progressBar) {
    window.addEventListener('scroll', onScroll, { passive: true });
    window.addEventListener('resize', onScroll, { passive: true });
    updateHeader(); // set initial state
  }

  // ── Scroll-reveal animations ──
  // Only enabled when the browser supports IntersectionObserver AND the user
  // has not requested reduced motion — otherwise content is left fully visible.
  var prefersReduced = window.matchMedia &&
    window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  if ('IntersectionObserver' in window && !prefersReduced) {
    var revealSelectors = [
      '.section-label', '.section-title', '.section-subtitle',
      '.lab-card', '.research-card', '.featured-tool-card',
      '.news-item', '.tool-card', '.pipe-form-section',
      '.analysis-card'
    ];
    var targets = document.querySelectorAll(revealSelectors.join(','));

    // Stagger items that share a direct parent so grids animate in sequence.
    var groupCounts = {};
    targets.forEach(function (el) {
      el.classList.add('js-reveal');
      var parent = el.parentElement;
      if (parent) {
        var key = parent.__revealKey || (parent.__revealKey = Math.random().toString(36).slice(2));
        var idx = groupCounts[key] = (groupCounts[key] || 0) + 1;
        if (idx <= 4) el.setAttribute('data-reveal-delay', idx);
      }
    });

    var observer = new IntersectionObserver(function (entries, obs) {
      entries.forEach(function (entry) {
        if (entry.isIntersecting) {
          entry.target.classList.add('is-visible');
          obs.unobserve(entry.target);
        }
      });
    }, { threshold: 0.12, rootMargin: '0px 0px -40px 0px' });

    targets.forEach(function (el) { observer.observe(el); });

    // Safety net: if anything is still hidden after 2s (e.g. observer quirk),
    // force it visible so content is never lost.
    setTimeout(function () {
      targets.forEach(function (el) { el.classList.add('is-visible'); });
    }, 2000);
  }

})();
