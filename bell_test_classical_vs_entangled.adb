-- -*- mode: indented-text; tab-width: 2; encoding: utf-8; -*-
----------------------------------------------------------------------

-- This is free and unencumbered software released into the public domain.
--
-- Anyone is free to copy, modify, publish, use, compile, sell, or
-- distribute this software, either in source code form or as a compiled
-- binary, for any purpose, commercial or non-commercial, and by any
-- means.
--
-- In jurisdictions that recognize copyright laws, the author or authors
-- of this software dedicate any and all copyright interest in the
-- software to the public domain. We make this dedication for the benefit
-- of the public at large and to the detriment of our heirs and
-- successors. We intend this dedication to be an overt act of
-- relinquishment in perpetuity of all present and future rights to this
-- software under copyright law.
--
-- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
-- EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
-- MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
-- IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
-- OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
-- ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
-- OTHER DEALINGS IN THE SOFTWARE.

----------------------------------------------------------------------

pragma ada_2022;
pragma wide_character_encoding (utf8);

with ada.assertions;
with ada.wide_wide_text_io;
with ada.numerics;
with ada.numerics.generic_elementary_functions;

procedure bell_test_classical_vs_entangled is

  subtype pair_range is integer range 1 .. 2;

  type scalar is digits 15;

  use ada.assertions;
  use ada.wide_wide_text_io;
  use ada.numerics;

  package scalar_elementary_functions is
    new ada.numerics.generic_elementary_functions (scalar);
  use scalar_elementary_functions;

  package scalar_io is new float_io (scalar);
  use scalar_io;

  π : constant scalar := pi;

----------------------------------------------------------------------

-- For the sake of reproducibility, let us write our own random number
-- generator. It will be a simple linear congruential generator. The
-- author has used one like it, in quicksorts and quickselects to
-- select the pivot. It is good enough for our purpose.

  type uint64 is mod 2 ** 64;

-- The multiplier lcg_a comes from Steele, Guy; Vigna, Sebastiano (28
-- September 2021). ‘Computationally easy, spectrally good multipliers
-- for congruential pseudorandom number generators’.
-- arXiv:2001.05304v3 [cs.DS]

  lcg_a : constant uint64 := 16#F1357AEA2E62A9C5#;

-- The value of lcg_c is not critical, but should be odd.

  lcg_c : constant uint64 := 1;

  seed  : uint64 := 0;

--
-- uniform_scalar: returns a non-negative scalar less than 1.
--
  function uniform_scalar
  return scalar
  with post => 0.0 <= uniform_scalar'result
                 and uniform_scalar'result < 1.0 is
    randval : scalar;
  begin
    -- Take the high 48 bits of the seed and divide it by 2**48.
    randval := scalar (seed / (2**16)) / scalar (2**48);

    -- Update the seed.
    seed := (lcg_a * seed) + lcg_c;

    return randval;
  end uniform_scalar;

----------------------------------------------------------------------

  type orientation is ('⇕', '⇔');
  type orientation_pair is array (pair_range) of orientation;

  subtype polarizing_beam_splitter is scalar;

  type plus_minus is ('⊕', '⊖');
  type plus_minus_pair is array (pair_range) of plus_minus;

  type event_record is
    record
      orientations : orientation_pair;
      detections   : plus_minus_pair;
    end record;

  function split_beam (φ : polarizing_beam_splitter;
                       σ : orientation)
  return plus_minus is
  begin
    return
      (case σ is
         when '⇕' =>
           (if uniform_scalar < sin (φ) ** 2 then '⊕' else '⊖'),
         when '⇔' =>
           (if uniform_scalar < cos (φ) ** 2 then '⊕' else '⊖'));
  end split_beam;

  function simulate_event_classically (φ1 : polarizing_beam_splitter;
                                       φ2 : polarizing_beam_splitter)
  return event_record is
    ev : event_record;
  begin
    ev.orientations := (if uniform_scalar < 0.5 then
                         ('⇕', '⇔')
                       else
                         ('⇔', '⇕'));
    ev.detections(1) := split_beam (φ1, ev.orientations(1));
    ev.detections(2) := split_beam (φ2, ev.orientations(2));
    return ev;
  end simulate_event_classically;

  function simulate_event_entangled (φ1 : polarizing_beam_splitter;
                                     φ2 : polarizing_beam_splitter)
  return event_record is
    ev : event_record;
  begin
    ev.orientations(1) := (if uniform_scalar < 0.5 then '⇕' else '⇔');
    ev.detections(1) := split_beam (φ1, ev.orientations(1));
    case ev.orientations(1) is
      when '⇕' =>
        ev.orientations(2) := '⇔';
        ev.detections(2) :=
          (if uniform_scalar < cos (φ2) ** 2 then '⊕' else '⊖');
      when '⇔' =>
        ev.orientations(2) := '⇕';
        ev.detections(2) :=
          (if uniform_scalar < sin (φ2) ** 2 then '⊕' else '⊖');
    end case;
    return ev;
  end simulate_event_entangled;

  type event_simulation is
    access function (φ1 : polarizing_beam_splitter;
                     φ2 : polarizing_beam_splitter)
           return event_record;

  type series_record is
    record
      φ1       : polarizing_beam_splitter;
      φ2       : polarizing_beam_splitter;
      num_ev   : positive;
      num_vhpp : natural;       -- vertical, horizontal, plus, minus
      num_vhpm : natural;
      num_vhmp : natural;
      num_vhmm : natural;
      num_hvpp : natural;
      num_hvpm : natural;
      num_hvmp : natural;
      num_hvmm : natural;
    end record;

  function simulate_event_series (ev_sim : event_simulation;
                                  φ1     : polarizing_beam_splitter;
                                  φ2     : polarizing_beam_splitter;
                                  num_ev : positive)
  return series_record is
    rec : series_record;
    ev  : event_record;
  begin
    rec.φ1 := φ1;
    rec.φ2 := φ2;
    rec.num_ev := num_ev;
    rec.num_vhpp := 0;
    rec.num_vhpm := 0;
    rec.num_vhmp := 0;
    rec.num_vhmm := 0;
    rec.num_hvpp := 0;
    rec.num_hvpm := 0;
    rec.num_hvmp := 0;
    rec.num_hvmm := 0;
    for i in 1 .. num_ev loop
      ev := ev_sim (φ1, φ2);
      case ev.orientations(1) is
        when '⇕' =>
          assert (ev.orientations(2) = '⇔');
          case ev.detections(1) is
            when '⊕' =>
              case ev.detections(2) is
                when '⊕' =>
                  rec.num_vhpp := @ + 1;
                when '⊖' =>
                  rec.num_vhpm := @ + 1;
              end case;
            when '⊖' =>
              case ev.detections(2) is
                when '⊕' =>
                  rec.num_vhmp := @ + 1;
                when '⊖' =>
                  rec.num_vhmm := @ + 1;
              end case;
          end case;
        when '⇔' =>
          assert (ev.orientations(2) = '⇕');
          case ev.detections(1) is
            when '⊕' =>
              case ev.detections(2) is
                when '⊕' =>
                  rec.num_hvpp := @ + 1;
                when '⊖' =>
                  rec.num_hvpm := @ + 1;
              end case;
            when '⊖' =>
              case ev.detections(2) is
                when '⊕' =>
                  rec.num_hvmp := @ + 1;
                when '⊖' =>
                  rec.num_hvmm := @ + 1;
              end case;
          end case;
      end case;
    end loop;
    return rec;
  end simulate_event_series;

begin
  null;
end bell_test_classical_vs_entangled;

----------------------------------------------------------------------
