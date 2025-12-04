"""Local time agent helper using a Hugging Face model (optional).

This module provides a minimal "agent" that can answer plain questions and
explicitly handle requests for the local time using a safe, local tool call.

Notes / assumptions:
- You must provide a valid Hugging Face model repo id if you want the agent to
  produce free-form answers via an LLM (for example, a SmolLM2 model ID).
  Replace the placeholder model name in examples with the actual repo id.
- This code does not automatically download or cache credentials. Set up HF
  auth (HUGGINGFACE_HUB_TOKEN) if the model is private.
"""

from __future__ import annotations

import datetime
import re
import typing as t
from zoneinfo import ZoneInfo

try:
    # optional import; the module works without transformers if you only need the tool
    from transformers import AutoModelForCausalLM, AutoTokenizer, pipeline
except Exception:  # pragma: no cover - transformers optional
    AutoModelForCausalLM = None
    AutoTokenizer = None
    pipeline = None


def get_local_time() -> str:
    """Return a short human-readable local time string (safe, no external calls).

    Uses the system local timezone. This is the tool the agent will call when a
    user asks for the current local time.
    """
    now = datetime.datetime.now().astimezone()
    tz = now.tzname() or str(now.tzinfo)
    iso = now.isoformat()
    return f"Local time: {iso} ({tz})"


_TIME_Q_RE = re.compile(r"\b(time|local time|what time|current time)\b", re.I)

# Small, editable mapping of common city names to IANA timezone names.
# This is intentionally minimal — extend as needed.
CITY_TZ_MAP = {
    "london": "Europe/London",
    "paris": "Europe/Paris",
    "new york": "America/New_York",
    "nyc": "America/New_York",
    "san francisco": "America/Los_Angeles",
    "los angeles": "America/Los_Angeles",
    "tokyo": "Asia/Tokyo",
    "seoul": "Asia/Seoul",
    "sydney": "Australia/Sydney",
    "berlin": "Europe/Berlin",
    "mumbai": "Asia/Kolkata",
}


def list_supported_cities() -> t.List[str]:
    return sorted(CITY_TZ_MAP.keys())


_IN_CITY_RE = re.compile(r"\bin\s+([A-Za-z0-9 _-]+)$", re.I)


def parse_city_from_prompt(text: str) -> t.Optional[str]:
    """Return a normalized city name if the prompt contains 'in <city>' or 'at <city>'."""
    m = re.search(r"\b(?:in|at)\s+([A-Za-z0-9 _-]+)\b", text, re.I)
    if not m:
        return None
    city = m.group(1).strip().lower()
    # normalize common abbreviations
    city = city.replace("’", "'")
    return city


def get_time_for_city(city: str) -> str:
    """Return the current time in the given city if known, else an error string."""
    tz_name = CITY_TZ_MAP.get(city.lower())
    if not tz_name:
        return f"Unknown city or timezone for '{city}'. Use --list-cities to see supported names."
    try:
        now = datetime.datetime.now(ZoneInfo(tz_name))
        iso = now.isoformat()
        return f"Local time in {city.title()}: {iso} ({tz_name})"
    except Exception as e:
        return f"Error computing time for {city}: {e}"


def is_time_query(text: str) -> bool:
    return bool(_TIME_Q_RE.search(text))


def load_model_pipeline(model_name: str, device: int = -1, use_auth_token: t.Optional[str] = None):
    """Load a Hugging Face text-generation pipeline for the given model name.

    model_name: HF repo id (e.g. 'owner/model').
    device: -1 for CPU, >=0 for CUDA device index.
    """
    if pipeline is None:
        raise RuntimeError("transformers is not installed. Install transformers to use an LLM.")

    tok = AutoTokenizer.from_pretrained(model_name, use_auth_token=use_auth_token)
    model = AutoModelForCausalLM.from_pretrained(model_name, use_auth_token=use_auth_token)
    gen = pipeline("text-generation", model=model, tokenizer=tok, device=device)
    return gen


def respond(
    prompt: str,
    model_name: t.Optional[str] = None,
    model_pipeline=None,
    device: int = -1,
    **gen_kwargs,
) -> str:
    """Produce a reply for `prompt`.

    Behavior:
    - If the prompt looks like a local-time question, return the direct tool output
      from `get_local_time()` (guaranteed accurate and local).
    - Otherwise, if `model_pipeline` or `model_name` is provided, use the model to
      generate a textual response. If neither is provided, return a short fallback.

    gen_kwargs are passed to the generation pipeline (max_new_tokens, etc.).
    """
    # If the prompt asks about time, check for explicit city mention
    if is_time_query(prompt):
        city = parse_city_from_prompt(prompt)
        if city:
            return get_time_for_city(city)
        return get_local_time()

    # If a pipeline instance is provided, use it directly
    if model_pipeline is not None:
        out = model_pipeline(prompt, **gen_kwargs)
        # pipeline returns list of dicts with 'generated_text'
        return out[0]["generated_text"]

    # If a model name is provided, lazily load a pipeline and run it
    if model_name:
        gen = load_model_pipeline(model_name, device=device)
        out = gen(prompt, **gen_kwargs)
        return out[0]["generated_text"]

    # Fallback: echo with a hint
    return f"(no LLM configured) I can answer time questions: {get_local_time()}"


if __name__ == "__main__":
    import argparse
    import readline

    p = argparse.ArgumentParser(description="Local-time aware chatbot (optional HF LLM)")
    p.add_argument("--model", help="Hugging Face model repo id to use for general answers")
    p.add_argument("--cpu", action="store_true", help="Force CPU usage for model (device=-1)")
    p.add_argument("--max-tokens", type=int, default=128)
    p.add_argument("--list-cities", action="store_true", help="Print supported cities and exit")
    args = p.parse_args()

    if args.list_cities:
        print("Supported cities:")
        for c in list_supported_cities():
            print(" -", c)
        raise SystemExit(0)

    device = -1 if args.cpu else -1

    # Simple REPL / chatbot loop
    history: t.List[t.Tuple[str, str]] = []  # (user, assistant)
    print("Local-time chatbot. Type 'exit' or Ctrl-C to quit.")
    while True:
        try:
            user = input("You: ").strip()
        except (KeyboardInterrupt, EOFError):
            print('\nGoodbye')
            break
        if not user:
            continue
        if user.lower() in ("exit", "quit"):
            print('Goodbye')
            break

        # Quick routing: if it's a time query and mentions a city, resolve via tool
        if is_time_query(user):
            city = parse_city_from_prompt(user)
            if city:
                out = get_time_for_city(city)
                print("Assistant:", out)
                history.append((user, out))
                continue

        # Otherwise, try to call the LLM if configured
        if args.model or 'model_pipeline' in globals():
            try:
                # naive conversation concatenation
                prompt = "\n".join([f"User: {u}\nAssistant: {a}" for u, a in history] + [f"User: {user}\nAssistant:"])
                resp = respond(prompt, model_name=args.model, device=device, max_new_tokens=args.max_tokens)
                print("Assistant:", resp)
                history.append((user, resp))
            except Exception as e:
                print("Assistant: Error generating response:", e)
        else:
            # no model configured — fallback
            out = respond(user)
            print("Assistant:", out)
            history.append((user, out))
