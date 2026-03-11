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

})();
