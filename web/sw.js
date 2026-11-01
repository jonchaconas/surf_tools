// Offline cache for the Predictive Surf Atlas. Data files are big and change only when re-exported,
// so they're served from the cache and refreshed in the background (stale-while-revalidate); the
// page itself and the forecast go to the network first and fall back to the cache when offline.
const CACHE = "psa-v1";
const FRESH = [/\/$/, /index\.html$/, /forecast\.json$/, /fcmap\.(json|b64\.txt)$/, /meta\.json$/];

self.addEventListener("install", e => { self.skipWaiting(); e.waitUntil(caches.open(CACHE).then(c => c.addAll(["./"]))); });
self.addEventListener("activate", e => {
  e.waitUntil(caches.keys().then(keys => Promise.all(keys.filter(k => k !== CACHE).map(k => caches.delete(k)))).then(() => self.clients.claim()));
});
self.addEventListener("fetch", e => {
  const url = new URL(e.request.url);
  if (e.request.method !== "GET" || url.origin !== location.origin) return;
  if (FRESH.some(r => r.test(url.pathname))) {
    e.respondWith(fetch(e.request).then(r => { const copy = r.clone(); caches.open(CACHE).then(c => c.put(e.request, copy)); return r; })
      .catch(() => caches.match(e.request)));
    return;
  }
  e.respondWith(caches.open(CACHE).then(async c => {
    const hit = await c.match(e.request);
    const net = fetch(e.request).then(r => { if (r.ok) c.put(e.request, r.clone()); return r; }).catch(() => hit);
    return hit || net;
  }));
});
