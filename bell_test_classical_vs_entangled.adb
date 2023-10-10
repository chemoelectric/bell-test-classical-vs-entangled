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

-- Simulation of both a classical model and an ‘entangled’ model of
-- a two-channel Bell test experiment with photons and polarizing
-- beam splitters.
--
-- Author: Barry Schwartz
-- Date completed: 9 October 2023
-- Revision history: See the public repository at
--   https://github.com/chemoelectric/bell-test-classical-vs-entangled

----------------------------------------------------------------------

pragma ada_2022;
pragma wide_character_encoding (utf8);
pragma assertion_policy (check);

with ada.assertions;
with ada.text_io;
with ada.wide_wide_text_io;
with ada.command_line;
with ada.numerics;
with ada.numerics.generic_elementary_functions;

procedure bell_test_classical_vs_entangled is

  subtype pair_range is integer range 1 .. 2;

  type scalar is digits 15;

  use ada.assertions;
  use ada.wide_wide_text_io;
  use ada.command_line;
  use ada.numerics;

  package tio renames ada.text_io;

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

  function simulate_event_classical (φ1 : polarizing_beam_splitter;
                                     φ2 : polarizing_beam_splitter)
  return event_record is

    --
    -- The ‘classical’ simulation goes as follows:
    --
    -- * A pair of perpendicularly polarized photons springs into
    --   existence.
    --
    -- * Each photon goes through a polarizing beam splitter and is
    --   directed towards its respective ⊕ photodetector or ⊖
    --   photodetector. The beam splitter obeys the Law of Malus (if
    --   we regard light intensity as a statistical phenomenon).
    --

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
    ev : event_record;
  begin
    ev.orientations := (if uniform_scalar < 0.5 then
                         ('⇕', '⇔')
                       else
                         ('⇔', '⇕'));
    ev.detections(1) := split_beam (φ1, ev.orientations(1));
    ev.detections(2) := split_beam (φ2, ev.orientations(2));
    return ev;
  end simulate_event_classical;

  type outcome_probabilities_record is
    record
      vhpp : scalar;            -- vertical, horizontal, plus, minus
      vhpm : scalar;
      vhmp : scalar;
      vhmm : scalar;
      hvpp : scalar;
      hvpm : scalar;
      hvmp : scalar;
      hvmm : scalar;
    end record;

  function outcome_probabilities (φ1 : polarizing_beam_splitter;
                                  φ2 : polarizing_beam_splitter)
  return outcome_probabilities_record is
    cos2φ1, cos2φ2   : scalar;
    sin2φ1, sin2φ2   : scalar;

    cos2φ1cos2φ2     : scalar;
    cos2φ1sin2φ2     : scalar;
    sin2φ1cos2φ2     : scalar;
    sin2φ1sin2φ2     : scalar;
  begin
    cos2φ1 := cos (φ1) ** 2;
    cos2φ2 := cos (φ2) ** 2;
    sin2φ1 := sin (φ1) ** 2;
    sin2φ2 := sin (φ2) ** 2;

    cos2φ1cos2φ2 := cos2φ1 * cos2φ2;
    cos2φ1sin2φ2 := cos2φ1 * sin2φ2;
    sin2φ1cos2φ2 := sin2φ1 * cos2φ2;
    sin2φ1sin2φ2 := sin2φ1 * sin2φ2;

    return (vhpp => 0.5 * sin2φ1cos2φ2,
            vhpm => 0.5 * sin2φ1sin2φ2,
            vhmp => 0.5 * cos2φ1cos2φ2,
            vhmm => 0.5 * cos2φ1sin2φ2,
            hvpp => 0.5 * cos2φ1sin2φ2,
            hvpm => 0.5 * cos2φ1cos2φ2,
            hvmp => 0.5 * sin2φ1sin2φ2,
            hvmm => 0.5 * sin2φ1cos2φ2);
  end outcome_probabilities;

  function simulate_event_entangled (φ1 : polarizing_beam_splitter;
                                     φ2 : polarizing_beam_splitter)
  return event_record is

    --
    -- The ‘entangled’ simulation goes as follows:
    --
    -- * The cumulative probability distribution of the cartesian
    --   product of (⇕, ⇔, ⊕, ⊖) is written down.
    --
    -- * An number is picked haphazardly from the interval [0,1].
    --
    -- * That number, along with the cumulative distribution,
    --   determines which of the possible outcomes is the one that
    --   will actually be.
    --
    -- (That this will produce the same statistics as the ‘classical’
    -- simulation ought to be obvious. ‘Entanglement’ is physicists
    -- imagining the probabilities of outcomes as the actual physical
    -- objects traveling through space, and in this regard is an
    -- example of what E. T. Jaynes called the ‘Mind Projection
    -- Fallacy’. Jaynes considered this fallacy an arrogance, but I
    -- think it ought to be considered a bona fide symptom of cerebral
    -- dysfunction, usually induced by a university education. Current
    -- education in quantum mechanics literally breaks a person’s
    -- brain, although hopefully not irreparably in too many cases.)
    --

    ev         : event_record;
    probs      : outcome_probabilities_record;
    cumulative : array (1 .. 8) of scalar;
    r          : scalar;
  begin
    probs := outcome_probabilities (φ1, φ2);

    cumulative(1) := probs.vhpp;
    cumulative(2) := cumulative(1) + probs.vhpm;
    cumulative(3) := cumulative(2) + probs.vhmp;
    cumulative(4) := cumulative(3) + probs.vhmm;
    cumulative(5) := cumulative(4) + probs.hvpp;
    cumulative(6) := cumulative(5) + probs.hvpm;
    cumulative(7) := cumulative(6) + probs.hvmp;
    cumulative(8) := cumulative(7) + probs.hvmm;

    assert (abs (cumulative(8) - 1.0) < 5.0e3 * scalar'model_epsilon);

    r := uniform_scalar;
    if r < cumulative(1) then
      ev.orientations := ('⇕', '⇔');
      ev.detections := ('⊕', '⊕');
    elsif r < cumulative(2) then
      ev.orientations := ('⇕', '⇔');
      ev.detections := ('⊕', '⊖');
    elsif r < cumulative(3) then
      ev.orientations := ('⇕', '⇔');
      ev.detections := ('⊖', '⊕');
    elsif r < cumulative(4) then
      ev.orientations := ('⇕', '⇔');
      ev.detections := ('⊖', '⊖');
    elsif r < cumulative(5) then
      ev.orientations := ('⇔', '⇕');
      ev.detections := ('⊕', '⊕');
    elsif r < cumulative(6) then
      ev.orientations := ('⇔', '⇕');
      ev.detections := ('⊕', '⊖');
    elsif r < cumulative(7) then
      ev.orientations := ('⇔', '⇕');
      ev.detections := ('⊖', '⊕');
    else
      ev.orientations := ('⇔', '⇕');
      ev.detections := ('⊖', '⊖');
    end if;

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

  function measure_correlation_coefficient (rec : series_record)
  return scalar is

    --
    -- This function computes an estimate of the correlation
    -- coefficient -cos(2(φ1-φ2)) from experimental data, for
    -- comparison with the nominal value and with other sets of
    -- experimental data.
    --
    -- Unlike ‘inequalities’, the correlation coefficient is an
    -- unimpeachable measure of correlation. If the measured values
    -- are reliably close to the nominal values, and the entire data
    -- set is employed (NO POST-SELECTION!) then the correlations ARE
    -- there.
    --

    freq_vhpp  : scalar;
    freq_vhpm  : scalar;
    freq_vhmp  : scalar;
    freq_vhmm  : scalar;
    freq_hvpp  : scalar;
    freq_hvpm  : scalar;
    freq_hvmp  : scalar;
    freq_hvmm  : scalar;
    cos2_cos2  : scalar;
    cos2_sin2  : scalar;
    sin2_cos2  : scalar;
    sin2_sin2  : scalar;
    cosφ1_sign : scalar;
    sinφ1_sign : scalar;
    cosφ2_sign : scalar;
    sinφ2_sign : scalar;
    cosφ1_φ2   : scalar;
    sinφ1_φ2   : scalar;
  begin
    -- Compute frequencies of events;
    freq_vhpp := scalar (rec.num_vhpp) / scalar (rec.num_ev);
    freq_vhpm := scalar (rec.num_vhpm) / scalar (rec.num_ev);
    freq_vhmp := scalar (rec.num_vhmp) / scalar (rec.num_ev);
    freq_vhmm := scalar (rec.num_vhmm) / scalar (rec.num_ev);
    freq_hvpp := scalar (rec.num_hvpp) / scalar (rec.num_ev);
    freq_hvpm := scalar (rec.num_hvpm) / scalar (rec.num_ev);
    freq_hvmp := scalar (rec.num_hvmp) / scalar (rec.num_ev);
    freq_hvmm := scalar (rec.num_hvmm) / scalar (rec.num_ev);

    -- Use these frequencies as estimates of products of squares of
    -- trigonometric functions (as one can infer from the experimental
    -- design).
    cos2_cos2 := freq_vhmp + freq_hvpm; -- cos²(φ1)×cos²(φ2)
    cos2_sin2 := freq_vhmm + freq_hvpp; -- cos²(φ1)×sin²(φ2)
    sin2_cos2 := freq_vhpp + freq_hvmm; -- sin²(φ1)×cos²(φ2)
    sin2_sin2 := freq_vhpm + freq_hvmp; -- sin²(φ1)×sin²(φ2)

    -- Because those quantities are squares, and we will need their
    -- square roots, there is a difficulty: there are two square roots
    -- for each square. To know which square roots to use, one needs
    -- to know which quadrants φ1 and φ2 are in. So let us figure out
    -- now which square roots to use. (This problem does not occur if
    -- all test angles are in Quadrant I, but in any case is not a
    -- major difficulty.)
    cosφ1_sign := (if cos (rec.φ1) < 0.0 then -1.0 else 1.0);
    sinφ1_sign := (if sin (rec.φ1) < 0.0 then -1.0 else 1.0);
    cosφ2_sign := (if cos (rec.φ2) < 0.0 then -1.0 else 1.0);
    sinφ2_sign := (if sin (rec.φ2) < 0.0 then -1.0 else 1.0);

    -- Use angle-difference identities to estimate cos(φ1-φ2) and
    -- sin(φ1-φ2).
    cosφ1_φ2 := (cosφ1_sign * cosφ2_sign * sqrt (cos2_cos2)) +
                   (sinφ1_sign * sinφ2_sign * sqrt (sin2_sin2));
    sinφ1_φ2 := (sinφ1_sign * cosφ2_sign * sqrt (sin2_cos2)) -
                   (cosφ1_sign * sinφ2_sign * sqrt (cos2_sin2));

    -- Return the estimate of -(cos²(φ1-φ2)-sin²(φ1-φ2)) =
    -- -cos(2(φ1-φ2)).
    return -((cosφ1_φ2 ** 2) - (sinφ1_φ2 ** 2));
  end measure_correlation_coefficient;

  φ1_deg, φ2_deg : scalar;
  φ1, φ2         : scalar;
  num_ev         : positive;
  rec_classical  : series_record;
  rec_entangled  : series_record;
  coef_classical : scalar;
  coef_entangled : scalar;
  probs          : outcome_probabilities_record;

  what_column      : constant positive_count := 4;
  nominal_column   : constant positive_count := 26;
  classical_column : constant positive_count := 39;
  entangled_column : constant positive_count := 52;

begin
  if argument_count /= 3 then
    put ("Usage: ");
    tio.put (command_name);
    put (" φ1 φ2 number_of_events");
    new_line;
    set_exit_status (1);
  else
  
    φ1_deg := scalar'value (argument(1));
    φ2_deg := scalar'value (argument(2));
    φ1 := φ1_deg * (π / 180.0);
    φ2 := φ2_deg * (π / 180.0);
    num_ev := positive'value (argument(3));
    rec_classical :=
      simulate_event_series (simulate_event_classical'access,
                             φ1, φ2, num_ev);
    rec_entangled :=
      simulate_event_series (simulate_event_entangled'access,
                             φ1, φ2, num_ev);
    coef_classical := measure_correlation_coefficient (rec_classical);
    coef_entangled := measure_correlation_coefficient (rec_entangled);
    probs := outcome_probabilities (φ1, φ2);

    new_line;

    set_col (to => what_column);
    put ("φ1 = "); put (φ1_deg, 5, 5, 0); put ("°"); new_line;
    set_col (to => what_column);
    put ("φ2 = "); put (φ2_deg, 5, 5, 0); put ("°"); new_line;

    new_line;

    set_col (to => what_column);
    put ("number of events ="); tio.put (num_ev'image); new_line;

    new_line;

    set_col (to => nominal_column); put (" nominal");
    set_col (to => classical_column); put (" classical");
    set_col (to => entangled_column); put (" entangled");
    new_line;

    set_col (to => what_column); put ("⇕ ⇔ ⊕ ⊕ frequency");
    set_col (to => nominal_column); put (probs.vhpp, 2, 5, 0);
    set_col (to => classical_column);
    put (scalar (rec_classical.num_vhpp) / scalar (rec_classical.num_ev),
         2, 5, 0);
    set_col (to => entangled_column);
    put (scalar (rec_entangled.num_vhpp) / scalar (rec_entangled.num_ev),
         2, 5, 0);
    new_line;

    set_col (to => what_column); put ("⇕ ⇔ ⊕ ⊖ frequency");
    set_col (to => nominal_column); put (probs.vhpm, 2, 5, 0);
    set_col (to => classical_column);
    put (scalar (rec_classical.num_vhpm) / scalar (rec_classical.num_ev),
         2, 5, 0);
    set_col (to => entangled_column);
    put (scalar (rec_entangled.num_vhpm) / scalar (rec_entangled.num_ev),
         2, 5, 0);
    new_line;

    set_col (to => what_column); put ("⇕ ⇔ ⊖ ⊕ frequency");
    set_col (to => nominal_column); put (probs.vhmp, 2, 5, 0);
    set_col (to => classical_column);
    put (scalar (rec_classical.num_vhmp) / scalar (rec_classical.num_ev),
         2, 5, 0);
    set_col (to => entangled_column);
    put (scalar (rec_entangled.num_vhmp) / scalar (rec_entangled.num_ev),
         2, 5, 0);
    new_line;

    set_col (to => what_column); put ("⇕ ⇔ ⊖ ⊖ frequency");
    set_col (to => nominal_column); put (probs.vhmm, 2, 5, 0);
    set_col (to => classical_column);
    put (scalar (rec_classical.num_vhmm) / scalar (rec_classical.num_ev),
         2, 5, 0);
    set_col (to => entangled_column);
    put (scalar (rec_entangled.num_vhmm) / scalar (rec_entangled.num_ev),
         2, 5, 0);
    new_line;

    set_col (to => what_column); put ("⇔ ⇕ ⊕ ⊕ frequency");
    set_col (to => nominal_column); put (probs.hvpp, 2, 5, 0);
    set_col (to => classical_column);
    put (scalar (rec_classical.num_hvpp) / scalar (rec_classical.num_ev),
         2, 5, 0);
    set_col (to => entangled_column);
    put (scalar (rec_entangled.num_hvpp) / scalar (rec_entangled.num_ev),
         2, 5, 0);
    new_line;

    set_col (to => what_column); put ("⇔ ⇕ ⊕ ⊖ frequency");
    set_col (to => nominal_column); put (probs.hvpm, 2, 5, 0);
    set_col (to => classical_column);
    put (scalar (rec_classical.num_hvpm) / scalar (rec_classical.num_ev),
         2, 5, 0);
    set_col (to => entangled_column);
    put (scalar (rec_entangled.num_hvpm) / scalar (rec_entangled.num_ev),
         2, 5, 0);
    new_line;

    set_col (to => what_column); put ("⇔ ⇕ ⊖ ⊕ frequency");
    set_col (to => nominal_column); put (probs.hvmp, 2, 5, 0);
    set_col (to => classical_column);
    put (scalar (rec_classical.num_hvmp) / scalar (rec_classical.num_ev),
         2, 5, 0);
    set_col (to => entangled_column);
    put (scalar (rec_entangled.num_hvmp) / scalar (rec_entangled.num_ev),
         2, 5, 0);
    new_line;

    set_col (to => what_column); put ("⇔ ⇕ ⊖ ⊖ frequency");
    set_col (to => nominal_column); put (probs.hvmm, 2, 5, 0);
    set_col (to => classical_column);
    put (scalar (rec_classical.num_hvmm) / scalar (rec_classical.num_ev),
         2, 5, 0);
    set_col (to => entangled_column);
    put (scalar (rec_entangled.num_hvmm) / scalar (rec_entangled.num_ev),
         2, 5, 0);
    new_line;

    set_col (to => what_column); put ("correlation coef.");
    set_col (to => nominal_column);
    put (-cos (2.0 * (φ1 - φ2)), 2, 5, 0);
    set_col (to => classical_column); put (coef_classical, 2, 5, 0);
    set_col (to => entangled_column); put (coef_entangled, 2, 5, 0);
    new_line;
  end if;
end bell_test_classical_vs_entangled;

----------------------------------------------------------------------
