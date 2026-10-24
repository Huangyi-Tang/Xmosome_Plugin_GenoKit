/* GenoKit browser adapter — main-thread Wasm (WKWebView blocks blob Workers under CSP). */
(() => {
  let cancelled = false;
  const post = (name, value) => window.webkit?.messageHandlers?.[name]?.postMessage(value);

  function notifyComplete(data, request) {
    try {
      if (typeof window.onGenoKitAnalysisComplete === 'function') {
        window.onGenoKitAnalysisComplete(data);
      } else if (typeof window.onPlinkitAnalysisComplete === 'function') {
        window.onPlinkitAnalysisComplete(data);
      }
    } catch (_) {
      /* UI callback errors must not block host completion. */
    }
    const artifacts = [...(data.outputs || []), ...(data.report ? [data.report] : [])];
    post('xmosomeExecutionComplete', {
      protocol_version: 1,
      run_id: request.run_id,
      status: data.status,
      result: data.status === 'success' ? data : null,
      artifacts,
      error: data.error || null,
    });
  }

  window.xmosomeBrowserCancel = () => {
    cancelled = true;
  };

  window.xmosomeBrowserExecute = async payload => {
    const request = payload.request || {};
    cancelled = false;
    try {
      let bytes;
      if (window.__pluginWasmBase64) {
        const raw = atob(window.__pluginWasmBase64);
        bytes = new Uint8Array(raw.length);
        for (let i = 0; i < raw.length; i++) bytes[i] = raw.charCodeAt(i);
      } else {
        const response = await fetch('../../Runtime/analysis.wasm');
        if (!response.ok) throw new Error('Browser runtime unavailable');
        bytes = new Uint8Array(await response.arrayBuffer());
      }

      if (cancelled) {
        notifyComplete(
          { status: 'failure', error: { code: 'cancelled', message: 'Cancelled.' } },
          request
        );
        return;
      }

      const { instance } = await WebAssembly.instantiate(bytes, {});
      const api = instance.exports;
      if (!api.pk_json_alloc || !api.pk_analyze_json || !api.pk_json_result) {
        throw new Error('Missing pk_json_alloc / pk_analyze_json / pk_json_result exports.');
      }

      const input = new TextEncoder().encode(
        JSON.stringify({ arguments: payload.arguments || {}, generated_at: Date.now() })
      );
      const ptr = api.pk_json_alloc(input.length);
      new Uint8Array(api.memory.buffer, ptr, input.length).set(input);
      const length = api.pk_analyze_json(ptr, input.length);
      const text = new TextDecoder().decode(
        new Uint8Array(api.memory.buffer, api.pk_json_result(), length)
      );

      let data;
      try {
        data = JSON.parse(text);
      } catch (parseError) {
        data = {
          status: 'failure',
          error: { code: 'parse-error', message: text || String(parseError) },
        };
      }
      notifyComplete(data, request);
    } catch (error) {
      notifyComplete(
        {
          status: 'failure',
          error: {
            code: 'execution-failed',
            message: (error && error.message) || String(error),
          },
        },
        request
      );
    }
  };
})();
