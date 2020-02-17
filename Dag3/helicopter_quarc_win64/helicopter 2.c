/*
 * helicopter.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "helicopter".
 *
 * Model version              : 1.190
 * Simulink Coder version : 8.9 (R2015b) 13-Aug-2015
 * C source code generated on : Mon Feb 03 20:45:10 2020
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "helicopter.h"
#include "helicopter_private.h"
#include "helicopter_dt.h"

/* Block signals (auto storage) */
B_helicopter_T helicopter_B;

/* Continuous states */
X_helicopter_T helicopter_X;

/* Block states (auto storage) */
DW_helicopter_T helicopter_DW;

/* Real-time model */
RT_MODEL_helicopter_T helicopter_M_;
RT_MODEL_helicopter_T *const helicopter_M = &helicopter_M_;

/*
 * Writes out MAT-file header.  Returns success or failure.
 * Returns:
 *      0 - success
 *      1 - failure
 */
int_T rt_WriteMat4FileHeader(FILE *fp, int32_T m, int32_T n, const char *name)
{
  typedef enum { ELITTLE_ENDIAN, EBIG_ENDIAN } ByteOrder;

  int16_T one = 1;
  ByteOrder byteOrder = (*((int8_T *)&one)==1) ? ELITTLE_ENDIAN : EBIG_ENDIAN;
  int32_T type = (byteOrder == ELITTLE_ENDIAN) ? 0: 1000;
  int32_T imagf = 0;
  int32_T name_len = (int32_T)strlen(name) + 1;
  if ((fwrite(&type, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&m, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&n, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&imagf, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(&name_len, sizeof(int32_T), 1, fp) == 0) ||
      (fwrite(name, sizeof(char), name_len, fp) == 0)) {
    return(1);
  } else {
    return(0);
  }
}                                      /* end rt_WriteMat4FileHeader */

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 4;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  helicopter_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Model output function */
void helicopter_output(void)
{
  /* local block i/o variables */
  real_T rtb_HILReadEncoderTimebase_o1;
  real_T rtb_HILReadEncoderTimebase_o2;
  real_T rtb_HILReadEncoderTimebase_o3;
  real_T rtb_TmpSignalConversionAtToFile[7];
  real_T *lastU;
  real_T rtb_Frontgain;
  if (rtmIsMajorTimeStep(helicopter_M)) {
    /* set solver stop time */
    if (!(helicopter_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&helicopter_M->solverInfo,
                            ((helicopter_M->Timing.clockTickH0 + 1) *
        helicopter_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&helicopter_M->solverInfo,
                            ((helicopter_M->Timing.clockTick0 + 1) *
        helicopter_M->Timing.stepSize0 + helicopter_M->Timing.clockTickH0 *
        helicopter_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(helicopter_M)) {
    helicopter_M->Timing.t[0] = rtsiGetT(&helicopter_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(helicopter_M)) {
    /* S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

    /* S-Function Block: helicopter/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_read_encoder(helicopter_DW.HILReadEncoderTimebase_Task,
        1, &helicopter_DW.HILReadEncoderTimebase_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
      } else {
        rtb_HILReadEncoderTimebase_o1 =
          helicopter_DW.HILReadEncoderTimebase_Buffer[0];
        rtb_HILReadEncoderTimebase_o2 =
          helicopter_DW.HILReadEncoderTimebase_Buffer[1];
        rtb_HILReadEncoderTimebase_o3 =
          helicopter_DW.HILReadEncoderTimebase_Buffer[2];
      }
    }
  }

  /* FromWorkspace: '<S6>/From' */
  {
    real_T *pDataValues = (real_T *) helicopter_DW.From_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helicopter_DW.From_PWORK.TimePtr;
    int_T currTimeIndex = helicopter_DW.From_IWORK.PrevIndex;
    real_T t = helicopter_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicopter_DW.From_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          helicopter_B.From = pDataValues[currTimeIndex];
        } else {
          helicopter_B.From = pDataValues[currTimeIndex + 1];
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;
        d1 = pDataValues[TimeIndex];
        d2 = pDataValues[TimeIndex + 1];
        helicopter_B.From = (real_T) rtInterpolate(d1, d2, f1, f2);
        pDataValues += 141;
      }
    }
  }

  if (rtmIsMajorTimeStep(helicopter_M)) {
    /* Gain: '<S4>/Travel: Count to rad' incorporates:
     *  Gain: '<S4>/Travel_gain'
     */
    helicopter_B.TravelCounttorad = helicopter_P.travel_gain *
      rtb_HILReadEncoderTimebase_o1 * helicopter_P.TravelCounttorad_Gain;

    /* Gain: '<S14>/Gain' */
    helicopter_B.Gain = helicopter_P.Gain_Gain * helicopter_B.TravelCounttorad;

    /* Gain: '<S4>/Pitch: Count to rad' */
    helicopter_B.PitchCounttorad = helicopter_P.PitchCounttorad_Gain *
      rtb_HILReadEncoderTimebase_o2;

    /* Gain: '<S11>/Gain' */
    helicopter_B.Gain_i = helicopter_P.Gain_Gain_a *
      helicopter_B.PitchCounttorad;
  }

  /* Gain: '<S15>/Gain' incorporates:
   *  TransferFcn: '<S4>/Travel: Transfer Fcn'
   */
  helicopter_B.Gain_d = (helicopter_P.TravelTransferFcn_C *
    helicopter_X.TravelTransferFcn_CSTATE + helicopter_P.TravelTransferFcn_D *
    helicopter_B.TravelCounttorad) * helicopter_P.Gain_Gain_l;

  /* Gain: '<S12>/Gain' incorporates:
   *  TransferFcn: '<S4>/Pitch: Transfer Fcn'
   */
  helicopter_B.Gain_b = (helicopter_P.PitchTransferFcn_C *
    helicopter_X.PitchTransferFcn_CSTATE + helicopter_P.PitchTransferFcn_D *
    helicopter_B.PitchCounttorad) * helicopter_P.Gain_Gain_ae;
  if (rtmIsMajorTimeStep(helicopter_M)) {
    /* Gain: '<S4>/Elevation: Count to rad' incorporates:
     *  Gain: '<S4>/Elevation_gain'
     */
    helicopter_B.ElevationCounttorad = helicopter_P.elevation_gain *
      rtb_HILReadEncoderTimebase_o3 * helicopter_P.ElevationCounttorad_Gain;

    /* Gain: '<S9>/Gain' */
    helicopter_B.Gain_e = helicopter_P.Gain_Gain_lv *
      helicopter_B.ElevationCounttorad;

    /* Sum: '<Root>/Sum' incorporates:
     *  Constant: '<Root>/elavation_offset [deg]'
     */
    helicopter_B.Sum = helicopter_B.Gain_e +
      helicopter_P.elavation_offsetdeg_Value;
  }

  /* Gain: '<S10>/Gain' incorporates:
   *  TransferFcn: '<S4>/Elevation: Transfer Fcn'
   */
  helicopter_B.Gain_dg = (helicopter_P.ElevationTransferFcn_C *
    helicopter_X.ElevationTransferFcn_CSTATE +
    helicopter_P.ElevationTransferFcn_D * helicopter_B.ElevationCounttorad) *
    helicopter_P.Gain_Gain_n;

  /* Gain: '<S3>/Gain1' */
  helicopter_B.Gain1[0] = helicopter_P.Gain1_Gain * helicopter_B.Gain;
  helicopter_B.Gain1[1] = helicopter_P.Gain1_Gain * helicopter_B.Gain_d;
  helicopter_B.Gain1[2] = helicopter_P.Gain1_Gain * helicopter_B.Gain_i;
  helicopter_B.Gain1[3] = helicopter_P.Gain1_Gain * helicopter_B.Gain_b;
  helicopter_B.Gain1[4] = helicopter_P.Gain1_Gain * helicopter_B.Sum;
  helicopter_B.Gain1[5] = helicopter_P.Gain1_Gain * helicopter_B.Gain_dg;

  /* Sum: '<S1>/Sum' incorporates:
   *  Constant: '<S1>/travel_offset'
   */
  helicopter_B.Sum_g = helicopter_P.travel_offset_Value + helicopter_B.Gain1[0];

  /* FromWorkspace: '<S6>/From1' */
  {
    real_T *pDataValues = (real_T *) helicopter_DW.From1_PWORK.DataPtr;
    real_T *pTimeValues = (real_T *) helicopter_DW.From1_PWORK.TimePtr;
    int_T currTimeIndex = helicopter_DW.From1_IWORK.PrevIndex;
    real_T t = helicopter_M->Timing.t[0];

    /* Get index */
    if (t <= pTimeValues[0]) {
      currTimeIndex = 0;
    } else if (t >= pTimeValues[140]) {
      currTimeIndex = 139;
    } else {
      if (t < pTimeValues[currTimeIndex]) {
        while (t < pTimeValues[currTimeIndex]) {
          currTimeIndex--;
        }
      } else {
        while (t >= pTimeValues[currTimeIndex + 1]) {
          currTimeIndex++;
        }
      }
    }

    helicopter_DW.From1_IWORK.PrevIndex = currTimeIndex;

    /* Post output */
    {
      real_T t1 = pTimeValues[currTimeIndex];
      real_T t2 = pTimeValues[currTimeIndex + 1];
      if (t1 == t2) {
        if (t < t1) {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter_B.From1[0])[elIdx] = pDataValues[currTimeIndex];
              pDataValues += 141;
            }
          }
        } else {
          {
            int_T elIdx;
            for (elIdx = 0; elIdx < 4; ++elIdx) {
              (&helicopter_B.From1[0])[elIdx] = pDataValues[currTimeIndex + 1];
              pDataValues += 141;
            }
          }
        }
      } else {
        real_T f1 = (t2 - t) / (t2 - t1);
        real_T f2 = 1.0 - f1;
        real_T d1;
        real_T d2;
        int_T TimeIndex= currTimeIndex;

        {
          int_T elIdx;
          for (elIdx = 0; elIdx < 4; ++elIdx) {
            d1 = pDataValues[TimeIndex];
            d2 = pDataValues[TimeIndex + 1];
            (&helicopter_B.From1[0])[elIdx] = (real_T) rtInterpolate(d1, d2, f1,
              f2);
            pDataValues += 141;
          }
        }
      }
    }
  }

  /* Gain: '<S6>/Feedback gain' incorporates:
   *  Sum: '<S6>/Subtract'
   */
  helicopter_B.Feedbackgain = (((helicopter_B.Sum_g - helicopter_B.From1[0]) *
    helicopter_P.K_lqr[0] + (helicopter_B.Gain1[1] - helicopter_B.From1[1]) *
    helicopter_P.K_lqr[1]) + (helicopter_B.Gain1[2] - helicopter_B.From1[2]) *
    helicopter_P.K_lqr[2]) + (helicopter_B.Gain1[3] - helicopter_B.From1[3]) *
    helicopter_P.K_lqr[3];

  /* Sum: '<S6>/Subtract1' */
  helicopter_B.Subtract1 = helicopter_B.From - helicopter_B.Feedbackgain;
  if (rtmIsMajorTimeStep(helicopter_M)) {
    /* SignalConversion: '<Root>/TmpSignal ConversionAtTo FileInport1' */
    rtb_TmpSignalConversionAtToFile[0] = helicopter_B.Subtract1;
    rtb_TmpSignalConversionAtToFile[1] = helicopter_B.Gain;
    rtb_TmpSignalConversionAtToFile[2] = helicopter_B.Gain_d;
    rtb_TmpSignalConversionAtToFile[3] = helicopter_B.Gain_i;
    rtb_TmpSignalConversionAtToFile[4] = helicopter_B.Gain_b;
    rtb_TmpSignalConversionAtToFile[5] = helicopter_B.Sum;
    rtb_TmpSignalConversionAtToFile[6] = helicopter_B.Gain_dg;

    /* ToFile: '<Root>/To File' */
    {
      if (!(++helicopter_DW.ToFile_IWORK.Decimation % 1) &&
          (helicopter_DW.ToFile_IWORK.Count*8)+1 < 100000000 ) {
        FILE *fp = (FILE *) helicopter_DW.ToFile_PWORK.FilePtr;
        if (fp != (NULL)) {
          real_T u[8];
          helicopter_DW.ToFile_IWORK.Decimation = 0;
          u[0] = helicopter_M->Timing.t[1];
          u[1] = rtb_TmpSignalConversionAtToFile[0];
          u[2] = rtb_TmpSignalConversionAtToFile[1];
          u[3] = rtb_TmpSignalConversionAtToFile[2];
          u[4] = rtb_TmpSignalConversionAtToFile[3];
          u[5] = rtb_TmpSignalConversionAtToFile[4];
          u[6] = rtb_TmpSignalConversionAtToFile[5];
          u[7] = rtb_TmpSignalConversionAtToFile[6];
          if (fwrite(u, sizeof(real_T), 8, fp) != 8) {
            rtmSetErrorStatus(helicopter_M,
                              "Error writing to MAT-file states.mat");
            return;
          }

          if (((++helicopter_DW.ToFile_IWORK.Count)*8)+1 >= 100000000) {
            (void)fprintf(stdout,
                          "*** The ToFile block will stop logging data before\n"
                          "    the simulation has ended, because it has reached\n"
                          "    the maximum number of elements (100000000)\n"
                          "    allowed in MAT-file states.mat.\n");
          }
        }
      }
    }
  }

  /* Sum: '<S5>/Sum1' incorporates:
   *  Constant: '<S5>/Vd_bias'
   *  Gain: '<S8>/K_pd'
   *  Gain: '<S8>/K_pp'
   *  Sum: '<S8>/Sum2'
   *  Sum: '<S8>/Sum3'
   */
  helicopter_B.Sum1 = ((helicopter_B.Subtract1 - helicopter_B.Gain1[2]) *
                       helicopter_P.K_pp - helicopter_P.K_pd *
                       helicopter_B.Gain1[3]) + helicopter_P.Vd_ff;
  if (rtmIsMajorTimeStep(helicopter_M)) {
  }

  /* Integrator: '<S7>/Integrator' */
  /* Limited  Integrator  */
  if (helicopter_X.Integrator_CSTATE >= helicopter_P.Integrator_UpperSat) {
    helicopter_X.Integrator_CSTATE = helicopter_P.Integrator_UpperSat;
  } else {
    if (helicopter_X.Integrator_CSTATE <= helicopter_P.Integrator_LowerSat) {
      helicopter_X.Integrator_CSTATE = helicopter_P.Integrator_LowerSat;
    }
  }

  /* Sum: '<S7>/Sum' incorporates:
   *  Constant: '<S5>/elevation_ref'
   */
  rtb_Frontgain = helicopter_P.elevation_ref_Value - helicopter_B.Gain1[4];

  /* Sum: '<S5>/Sum2' incorporates:
   *  Constant: '<S5>/Vs_bias'
   *  Gain: '<S7>/K_ed'
   *  Gain: '<S7>/K_ep'
   *  Integrator: '<S7>/Integrator'
   *  Sum: '<S7>/Sum1'
   */
  helicopter_B.Sum2 = ((helicopter_P.K_ep * rtb_Frontgain +
                        helicopter_X.Integrator_CSTATE) - helicopter_P.K_ed *
                       helicopter_B.Gain1[5]) + helicopter_P.Vs_ff;
  if (rtmIsMajorTimeStep(helicopter_M)) {
  }

  /* Gain: '<S7>/K_ei' */
  helicopter_B.K_ei = helicopter_P.K_ei * rtb_Frontgain;
  if (rtmIsMajorTimeStep(helicopter_M)) {
  }

  /* Derivative: '<S4>/Derivative' */
  if ((helicopter_DW.TimeStampA >= helicopter_M->Timing.t[0]) &&
      (helicopter_DW.TimeStampB >= helicopter_M->Timing.t[0])) {
    rtb_Frontgain = 0.0;
  } else {
    rtb_Frontgain = helicopter_DW.TimeStampA;
    lastU = &helicopter_DW.LastUAtTimeA;
    if (helicopter_DW.TimeStampA < helicopter_DW.TimeStampB) {
      if (helicopter_DW.TimeStampB < helicopter_M->Timing.t[0]) {
        rtb_Frontgain = helicopter_DW.TimeStampB;
        lastU = &helicopter_DW.LastUAtTimeB;
      }
    } else {
      if (helicopter_DW.TimeStampA >= helicopter_M->Timing.t[0]) {
        rtb_Frontgain = helicopter_DW.TimeStampB;
        lastU = &helicopter_DW.LastUAtTimeB;
      }
    }

    rtb_Frontgain = (helicopter_B.PitchCounttorad - *lastU) /
      (helicopter_M->Timing.t[0] - rtb_Frontgain);
  }

  /* End of Derivative: '<S4>/Derivative' */

  /* Gain: '<S13>/Gain' */
  helicopter_B.Gain_l = helicopter_P.Gain_Gain_a1 * rtb_Frontgain;
  if (rtmIsMajorTimeStep(helicopter_M)) {
  }

  /* Gain: '<S2>/Back gain' incorporates:
   *  Sum: '<S2>/Subtract'
   */
  rtb_Frontgain = (helicopter_B.Sum2 - helicopter_B.Sum1) *
    helicopter_P.Backgain_Gain;

  /* Saturate: '<S4>/Back motor: Saturation' */
  if (rtb_Frontgain > helicopter_P.BackmotorSaturation_UpperSat) {
    helicopter_B.BackmotorSaturation = helicopter_P.BackmotorSaturation_UpperSat;
  } else if (rtb_Frontgain < helicopter_P.BackmotorSaturation_LowerSat) {
    helicopter_B.BackmotorSaturation = helicopter_P.BackmotorSaturation_LowerSat;
  } else {
    helicopter_B.BackmotorSaturation = rtb_Frontgain;
  }

  /* End of Saturate: '<S4>/Back motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_M)) {
  }

  /* Gain: '<S2>/Front gain' incorporates:
   *  Sum: '<S2>/Add'
   */
  rtb_Frontgain = (helicopter_B.Sum1 + helicopter_B.Sum2) *
    helicopter_P.Frontgain_Gain;

  /* Saturate: '<S4>/Front motor: Saturation' */
  if (rtb_Frontgain > helicopter_P.FrontmotorSaturation_UpperSat) {
    helicopter_B.FrontmotorSaturation =
      helicopter_P.FrontmotorSaturation_UpperSat;
  } else if (rtb_Frontgain < helicopter_P.FrontmotorSaturation_LowerSat) {
    helicopter_B.FrontmotorSaturation =
      helicopter_P.FrontmotorSaturation_LowerSat;
  } else {
    helicopter_B.FrontmotorSaturation = rtb_Frontgain;
  }

  /* End of Saturate: '<S4>/Front motor: Saturation' */
  if (rtmIsMajorTimeStep(helicopter_M)) {
    /* S-Function (hil_write_analog_block): '<S4>/HIL Write Analog' */

    /* S-Function Block: helicopter/Helicopter_interface/HIL Write Analog (hil_write_analog_block) */
    {
      t_error result;
      helicopter_DW.HILWriteAnalog_Buffer[0] = helicopter_B.FrontmotorSaturation;
      helicopter_DW.HILWriteAnalog_Buffer[1] = helicopter_B.BackmotorSaturation;
      result = hil_write_analog(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILWriteAnalog_channels, 2,
        &helicopter_DW.HILWriteAnalog_Buffer[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
      }
    }
  }
}

/* Model update function */
void helicopter_update(void)
{
  real_T *lastU;

  /* Update for Derivative: '<S4>/Derivative' */
  if (helicopter_DW.TimeStampA == (rtInf)) {
    helicopter_DW.TimeStampA = helicopter_M->Timing.t[0];
    lastU = &helicopter_DW.LastUAtTimeA;
  } else if (helicopter_DW.TimeStampB == (rtInf)) {
    helicopter_DW.TimeStampB = helicopter_M->Timing.t[0];
    lastU = &helicopter_DW.LastUAtTimeB;
  } else if (helicopter_DW.TimeStampA < helicopter_DW.TimeStampB) {
    helicopter_DW.TimeStampA = helicopter_M->Timing.t[0];
    lastU = &helicopter_DW.LastUAtTimeA;
  } else {
    helicopter_DW.TimeStampB = helicopter_M->Timing.t[0];
    lastU = &helicopter_DW.LastUAtTimeB;
  }

  *lastU = helicopter_B.PitchCounttorad;

  /* End of Update for Derivative: '<S4>/Derivative' */
  if (rtmIsMajorTimeStep(helicopter_M)) {
    rt_ertODEUpdateContinuousStates(&helicopter_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++helicopter_M->Timing.clockTick0)) {
    ++helicopter_M->Timing.clockTickH0;
  }

  helicopter_M->Timing.t[0] = rtsiGetSolverStopTime(&helicopter_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.002s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++helicopter_M->Timing.clockTick1)) {
      ++helicopter_M->Timing.clockTickH1;
    }

    helicopter_M->Timing.t[1] = helicopter_M->Timing.clockTick1 *
      helicopter_M->Timing.stepSize1 + helicopter_M->Timing.clockTickH1 *
      helicopter_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Derivatives for root system: '<Root>' */
void helicopter_derivatives(void)
{
  boolean_T lsat;
  boolean_T usat;
  XDot_helicopter_T *_rtXdot;
  _rtXdot = ((XDot_helicopter_T *) helicopter_M->ModelData.derivs);

  /* Derivatives for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  _rtXdot->TravelTransferFcn_CSTATE = 0.0;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter_P.TravelTransferFcn_A *
    helicopter_X.TravelTransferFcn_CSTATE;
  _rtXdot->TravelTransferFcn_CSTATE += helicopter_B.TravelCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  _rtXdot->PitchTransferFcn_CSTATE = 0.0;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_P.PitchTransferFcn_A *
    helicopter_X.PitchTransferFcn_CSTATE;
  _rtXdot->PitchTransferFcn_CSTATE += helicopter_B.PitchCounttorad;

  /* Derivatives for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  _rtXdot->ElevationTransferFcn_CSTATE = 0.0;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter_P.ElevationTransferFcn_A *
    helicopter_X.ElevationTransferFcn_CSTATE;
  _rtXdot->ElevationTransferFcn_CSTATE += helicopter_B.ElevationCounttorad;

  /* Derivatives for Integrator: '<S7>/Integrator' */
  lsat = (helicopter_X.Integrator_CSTATE <= helicopter_P.Integrator_LowerSat);
  usat = (helicopter_X.Integrator_CSTATE >= helicopter_P.Integrator_UpperSat);
  if (((!lsat) && (!usat)) || (lsat && (helicopter_B.K_ei > 0.0)) || (usat &&
       (helicopter_B.K_ei < 0.0))) {
    _rtXdot->Integrator_CSTATE = helicopter_B.K_ei;
  } else {
    /* in saturation */
    _rtXdot->Integrator_CSTATE = 0.0;
  }

  /* End of Derivatives for Integrator: '<S7>/Integrator' */
}

/* Model initialize function */
void helicopter_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q8_usb", "0", &helicopter_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options(helicopter_DW.HILInitialize_Card,
      "update_rate=normal;decimation=1", 32);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear(helicopter_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_M, _rt_error_message);
      return;
    }

    if ((helicopter_P.HILInitialize_set_analog_input_ && !is_switching) ||
        (helicopter_P.HILInitialize_set_analog_inpu_m && is_switching)) {
      {
        int_T i1;
        real_T *dw_AIMinimums = &helicopter_DW.HILInitialize_AIMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMinimums[i1] = helicopter_P.HILInitialize_analog_input_mini;
        }
      }

      {
        int_T i1;
        real_T *dw_AIMaximums = &helicopter_DW.HILInitialize_AIMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AIMaximums[i1] = helicopter_P.HILInitialize_analog_input_maxi;
        }
      }

      result = hil_set_analog_input_ranges(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_analog_input_chan, 8U,
        &helicopter_DW.HILInitialize_AIMinimums[0],
        &helicopter_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_P.HILInitialize_set_analog_output && !is_switching) ||
        (helicopter_P.HILInitialize_set_analog_outp_b && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOMinimums = &helicopter_DW.HILInitialize_AOMinimums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMinimums[i1] = helicopter_P.HILInitialize_analog_output_min;
        }
      }

      {
        int_T i1;
        real_T *dw_AOMaximums = &helicopter_DW.HILInitialize_AOMaximums[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOMaximums[i1] = helicopter_P.HILInitialize_analog_output_max;
        }
      }

      result = hil_set_analog_output_ranges(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_analog_output_cha, 8U,
        &helicopter_DW.HILInitialize_AOMinimums[0],
        &helicopter_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_P.HILInitialize_set_analog_outp_e && !is_switching) ||
        (helicopter_P.HILInitialize_set_analog_outp_j && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_P.HILInitialize_initial_analog_ou;
        }
      }

      result = hil_write_analog(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_analog_output_cha, 8U,
        &helicopter_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_P.HILInitialize_set_analog_outp_p) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_P.HILInitialize_watchdog_analog_o;
        }
      }

      result = hil_watchdog_set_analog_expiration_state
        (helicopter_DW.HILInitialize_Card,
         helicopter_P.HILInitialize_analog_output_cha, 8U,
         &helicopter_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_P.HILInitialize_set_encoder_param && !is_switching) ||
        (helicopter_P.HILInitialize_set_encoder_par_m && is_switching)) {
      {
        int_T i1;
        int32_T *dw_QuadratureModes =
          &helicopter_DW.HILInitialize_QuadratureModes[0];
        for (i1=0; i1 < 8; i1++) {
          dw_QuadratureModes[i1] = helicopter_P.HILInitialize_quadrature;
        }
      }

      result = hil_set_encoder_quadrature_mode(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_encoder_channels, 8U,
        (t_encoder_quadrature_mode *)
        &helicopter_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_P.HILInitialize_set_encoder_count && !is_switching) ||
        (helicopter_P.HILInitialize_set_encoder_cou_k && is_switching)) {
      {
        int_T i1;
        int32_T *dw_InitialEICounts =
          &helicopter_DW.HILInitialize_InitialEICounts[0];
        for (i1=0; i1 < 8; i1++) {
          dw_InitialEICounts[i1] = helicopter_P.HILInitialize_initial_encoder_c;
        }
      }

      result = hil_set_encoder_counts(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_encoder_channels, 8U,
        &helicopter_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_P.HILInitialize_set_pwm_params_at && !is_switching) ||
        (helicopter_P.HILInitialize_set_pwm_params__f && is_switching)) {
      uint32_T num_duty_cycle_modes = 0;
      uint32_T num_frequency_modes = 0;

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter_P.HILInitialize_pwm_modes;
        }
      }

      result = hil_set_pwm_mode(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_pwm_channels, 8U, (t_pwm_mode *)
        &helicopter_DW.HILInitialize_POModeValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        const uint32_T *p_HILInitialize_pwm_channels =
          helicopter_P.HILInitialize_pwm_channels;
        int32_T *dw_POModeValues = &helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          if (dw_POModeValues[i1] == PWM_DUTY_CYCLE_MODE || dw_POModeValues[i1] ==
              PWM_ONE_SHOT_MODE || dw_POModeValues[i1] == PWM_TIME_MODE) {
            helicopter_DW.HILInitialize_POSortedChans[num_duty_cycle_modes] =
              p_HILInitialize_pwm_channels[i1];
            helicopter_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes] =
              helicopter_P.HILInitialize_pwm_frequency;
            num_duty_cycle_modes++;
          } else {
            helicopter_DW.HILInitialize_POSortedChans[7U - num_frequency_modes] =
              p_HILInitialize_pwm_channels[i1];
            helicopter_DW.HILInitialize_POSortedFreqs[7U - num_frequency_modes] =
              helicopter_P.HILInitialize_pwm_frequency;
            num_frequency_modes++;
          }
        }
      }

      if (num_duty_cycle_modes > 0) {
        result = hil_set_pwm_frequency(helicopter_DW.HILInitialize_Card,
          &helicopter_DW.HILInitialize_POSortedChans[0], num_duty_cycle_modes,
          &helicopter_DW.HILInitialize_POSortedFreqs[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_M, _rt_error_message);
          return;
        }
      }

      if (num_frequency_modes > 0) {
        result = hil_set_pwm_duty_cycle(helicopter_DW.HILInitialize_Card,
          &helicopter_DW.HILInitialize_POSortedChans[num_duty_cycle_modes],
          num_frequency_modes,
          &helicopter_DW.HILInitialize_POSortedFreqs[num_duty_cycle_modes]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_M, _rt_error_message);
          return;
        }
      }

      {
        int_T i1;
        int32_T *dw_POModeValues = &helicopter_DW.HILInitialize_POModeValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POModeValues[i1] = helicopter_P.HILInitialize_pwm_configuration;
        }
      }

      {
        int_T i1;
        int32_T *dw_POAlignValues = &helicopter_DW.HILInitialize_POAlignValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POAlignValues[i1] = helicopter_P.HILInitialize_pwm_alignment;
        }
      }

      {
        int_T i1;
        int32_T *dw_POPolarityVals =
          &helicopter_DW.HILInitialize_POPolarityVals[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POPolarityVals[i1] = helicopter_P.HILInitialize_pwm_polarity;
        }
      }

      result = hil_set_pwm_configuration(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_pwm_channels, 8U,
        (t_pwm_configuration *) &helicopter_DW.HILInitialize_POModeValues[0],
        (t_pwm_alignment *) &helicopter_DW.HILInitialize_POAlignValues[0],
        (t_pwm_polarity *) &helicopter_DW.HILInitialize_POPolarityVals[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }

      {
        int_T i1;
        real_T *dw_POSortedFreqs = &helicopter_DW.HILInitialize_POSortedFreqs[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POSortedFreqs[i1] = helicopter_P.HILInitialize_pwm_leading_deadb;
        }
      }

      {
        int_T i1;
        real_T *dw_POValues = &helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_P.HILInitialize_pwm_trailing_dead;
        }
      }

      result = hil_set_pwm_deadband(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_pwm_channels, 8U,
        &helicopter_DW.HILInitialize_POSortedFreqs[0],
        &helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if ((helicopter_P.HILInitialize_set_pwm_outputs_a && !is_switching) ||
        (helicopter_P.HILInitialize_set_pwm_outputs_g && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_P.HILInitialize_initial_pwm_outpu;
        }
      }

      result = hil_write_pwm(helicopter_DW.HILInitialize_Card,
        helicopter_P.HILInitialize_pwm_channels, 8U,
        &helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }

    if (helicopter_P.HILInitialize_set_pwm_outputs_o) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_P.HILInitialize_watchdog_pwm_outp;
        }
      }

      result = hil_watchdog_set_pwm_expiration_state
        (helicopter_DW.HILInitialize_Card,
         helicopter_P.HILInitialize_pwm_channels, 8U,
         &helicopter_DW.HILInitialize_POValues[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(helicopter_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S4>/HIL Read Encoder Timebase' */

  /* S-Function Block: helicopter/Helicopter_interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader(helicopter_DW.HILInitialize_Card,
      helicopter_P.HILReadEncoderTimebase_samples_,
      helicopter_P.HILReadEncoderTimebase_channels, 3,
      &helicopter_DW.HILReadEncoderTimebase_Task);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(helicopter_M, _rt_error_message);
    }
  }

  /* Start for FromWorkspace: '<S6>/From' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.308996938995715,
      1.3089969389957086, 1.308996938995701, 1.3089969389956908,
      1.3089969389956766, 1.3089969389956564, 1.3089969389956264,
      1.3089969389955787, 1.3089969389954947, 1.3089969389953195,
      1.3089969389946767, 1.281129019932987, -0.1163319635607416,
      -1.1273342905627004, -1.3089969389922591, -1.3089969389938949,
      -1.3089969389941987, -1.3089969389942275, -1.3089969389940854,
      -1.3089969389937071, -1.3089969389928666, -1.3089969389912126,
      -1.3089969303676823, -1.1062106034146504, -0.85872814734843184,
      -0.63872737187212092, -0.44991582394205615, -0.29325225079622552,
      -0.16766299409831406, -0.070687812188273622, 0.00097366719177557307,
      0.051025248722447084, 0.0832443476664491, 0.10124805192083643,
      0.10834012484106399, 0.10742298167817117, 0.10095946118666388,
      0.0909709647507475, 0.079060670615851661, 0.06645269764036453,
      0.054040099935823219, 0.042436333751664668, 0.032026325260893641,
      0.023014489444279949, 0.015468028320421814, 0.0093545994704743868,
      0.00457402250151822, 0.00098411023686361364, -0.0015790002741695416,
      -0.0032854597715907545, -0.0043006065694800441, -0.0047774649130612138,
      -0.0048521657098171134, -0.0046416587825661307, -0.0042431514750470752,
      -0.0037347848980459988, -0.003177139016666445, -0.0026152358071331889,
      -0.0020807822026513234, -0.0015944592389942045, -0.0011681195383548728,
      -0.0008068017068214113, -0.000510507642220764, -0.00027571783773246627,
      -9.6641452663394993E-5, 3.3786742752603093E-5, 0.00012314037661078893,
      0.00017893843125710857, 0.00020825457578199899, 0.0002174681688968349,
      0.00021212779357903836, 0.00019690085849002561, 0.00017558610315826688,
      0.00015116940632081926, 0.00012590684828674465, 0.00010142232959997351,
      7.8810078212270061E-5, 5.87350182800818E-5, 4.1526199312069279E-5,
      2.7260298712230993E-5, 1.5833638709886872E-5, 7.0222385059341949E-6,
      5.3019990681332961E-7, -3.9727513695950385E-6, -6.8234269651606111E-6,
      -8.3457348491293E-6, -8.8382784386865775E-6, -8.567652359181375E-6,
      -7.76613115237074E-6, -6.6324771437316063E-6, -5.33458874845779E-6,
      -4.0126673600417535E-6, -2.7815351804703228E-6, -1.7307862758002741E-6,
      -9.2178643871066577E-7, -3.8144947066167149E-7, -9.4600184104393164E-8,
      -4.8223842571264433E-19, -4.0776667315588432E-19, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0 } ;

    helicopter_DW.From_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter_DW.From_PWORK.DataPtr = (void *) pDataValues0;
    helicopter_DW.From_IWORK.PrevIndex = 0;
  }

  /* Start for FromWorkspace: '<S6>/From1' */
  {
    static real_T pTimeValues0[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
      2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
      5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5, 7.75, 8.0, 8.25,
      8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 10.75, 11.0, 11.25,
      11.5, 11.75, 12.0, 12.25, 12.5, 12.75, 13.0, 13.25, 13.5, 13.75, 14.0,
      14.25, 14.5, 14.75, 15.0, 15.25, 15.5, 15.75, 16.0, 16.25, 16.5, 16.75,
      17.0, 17.25, 17.5, 17.75, 18.0, 18.25, 18.5, 18.75, 19.0, 19.25, 19.5,
      19.75, 20.0, 20.25, 20.5, 20.75, 21.0, 21.25, 21.5, 21.75, 22.0, 22.25,
      22.5, 22.75, 23.0, 23.25, 23.5, 23.75, 24.0, 24.25, 24.5, 24.75, 25.0,
      25.25, 25.5, 25.75, 26.0, 26.25, 26.5, 26.75, 27.0, 27.25, 27.5, 27.75,
      28.0, 28.25, 28.5, 28.75, 29.0, 29.25, 29.5, 29.75, 30.0, 30.25, 30.5,
      30.75, 31.0, 31.25, 31.5, 31.75, 32.0, 32.25, 32.5, 32.75, 33.0, 33.25,
      33.5, 33.75, 34.0, 34.25, 34.5, 34.75, 35.0 } ;

    static real_T pDataValues0[] = { 6.2831853071795862, 6.2831853071795862,
      6.2831853071795862, 6.2831853071795862, 6.2831853071795862,
      6.2831853071795862, 6.2831853071795862, 6.2831853071795862,
      6.2831853071795862, 6.2831853071795862, 6.2831853071795862,
      6.2831853071795862, 6.2831853071795862, 6.2831853071795862,
      6.2831853071795862, 6.2831853071795862, 6.2831853071795862,
      6.2831853071795862, 6.2831853071795862, 6.2831853071795862,
      6.2831853071795862, 6.2831853071795862, 6.2831853071795862,
      6.2831853071795862, 6.2738090266114188, 6.2438986915989645,
      6.1840677076654309, 6.0877765035420248, 5.95113842612675,
      5.7722422803047833, 5.550426704926811, 5.2856978261978016,
      4.9783373991786366, 4.6286781700167, 4.23699917827712, 3.8036940079221764,
      3.33912516660695, 2.8627024609439329, 2.394524360723703,
      1.9504529752528792, 1.5410233496432477, 1.1721894122155887,
      0.84666429395023768, 0.56520574046066818, 0.32758951967006783,
      0.13323433097742668, -0.018458517082349836, -0.12941277607019347,
      -0.20344195802275125, -0.24586164007919895, -0.26276719935335491,
      -0.26032834355214868, -0.24425703300142507, -0.21948470023976038,
      -0.190022594098239, -0.1589560731943622, -0.12852282983607236,
      -0.10023425368993027, -0.075011056320185621, -0.053315088325392744,
      -0.03526755873990188, -0.020749491208440836, -0.009483689269017671,
      -0.0010993540301828758, 0.0048186394828409525, 0.00869352565544429,
      0.010931242919535369, 0.011903767075009472, 0.011939339803649364,
      0.011317988939409456, 0.010270970344002955, 0.0089829813064224891,
      0.0075961998522417074, 0.0062153917422489817, 0.0049134957976574309,
      0.0037372471377995326, 0.002712526022003487, 0.0018492270336005116,
      0.0011455299309277749, 0.00059152096958471081, 0.00017216380405654522,
      -0.00013034549789301336, -0.000334781797068167, -0.00045958205337109681,
      -0.00052194853211780527, -0.00053728869249687866, -0.0005189208106202984,
      -0.00047798132729743816, -0.00042347823782655237, -0.00036244368441308751,
      -0.00030014763390883158, -0.00024034268656791986, -0.00018551739054784732,
      -0.00013714178556401311, -9.589421854714803E-5, -6.18627860772542E-5,
      -3.4718132566206136E-5, -1.3856869187213526E-5, 1.4833086528203187E-6,
      1.213451845734702E-5, 1.8932594091723334E-5, 2.2670645756809856E-5,
      2.4068869459744772E-5, 2.3757343658654898E-5, 2.2268828978252963E-5,
      2.0038943627269771E-5, 1.7411479767404522E-5, 1.4647022932548654E-5,
      1.1933415409740841E-5, 9.3969507403474466E-6, 7.11349175493342E-6,
      5.1189649209326488E-6, 3.4188984030006639E-6, 1.9968409252937627E-6,
      8.2162457060898882E-7, -1.464820608125009E-7, -9.5064318817800092E-7,
      -1.6332621640687294E-6, -2.2326806530518027E-6, -2.7807982517158371E-6,
      -3.3017183765099408E-6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, -0.037505122272670145, -0.11964134004981758,
      -0.23932393573413469, -0.38516481649362527, -0.54655230966109658,
      -0.71558458328786712, -0.88726230151189134, -1.0589155149160363,
      -1.2294417080766631, -1.398636916647747, -1.5667159669583175,
      -1.733220681419777, -1.8582753652609043, -1.9056908226520701,
      -1.8727124008809184, -1.7762855418832957, -1.6377185024385263,
      -1.4753357497106365, -1.3021004730614034, -1.1258342139582778,
      -0.95046488316240152, -0.7774207547705646, -0.606771392239106,
      -0.44381703595137456, -0.29611672781023118, -0.16967872822579069,
      -0.0676222370966239, 0.0097554232048250276, 0.064285242202894374,
      0.099089331046658707, 0.11784842456608552, 0.1242660836155072,
      0.12173297343315938, 0.11315430458456836, 0.10089278947897863,
      0.086783871979171492, 0.072190118341963469, 0.058072270125844171,
      0.045063207757692653, 0.033537340955339182, 0.023671974052095313,
      0.015499544690413349, 0.0089508690563643144, 0.0038900966218964081,
      0.00014229091455957458, -0.0024854034569596357, -0.0041880743816260059,
      -0.0051519561503218627, -0.0055471258167231283, -0.0055232324399709029,
      -0.0052075837783662057, -0.0047049946394315924, -0.0040988844631841825,
      -0.0034531959536119009, -0.0028147884106909469, -0.002216035845372256,
      -0.0016774286621126626, -0.0012100372077982343, -0.0008177451967006146,
      -0.00049920102521171911, -0.0002494659149868339, -6.1360641516293722E-5,
      7.3471527506320785E-5, 0.0001637579332914409, 0.00021801235788354326,
      0.00024413821365385944, 0.00024918420201702361, 0.00023921978936364678,
      0.00021930118408029016, 0.00019350241993533685, 0.00016499026806746038,
      0.00013612572987957525, 0.00010857861404419229, 8.3445053515970448E-5,
      6.1360711360135384E-5, 4.26048392181068E-5, 2.7192302537505254E-5,
      1.4952206660346095E-5, 5.5928948117396709E-6, -1.2461032043594969E-6,
      -5.9540587216077368E-6, -8.919541403932769E-6, -1.0509855439460997E-5,
      -1.1057827339423475E-5, -1.085443009123126E-5, -1.0145858677573575E-5,
      -9.133835941656107E-6, -7.978107336003086E-6, -6.8002660717279395E-6,
      -5.6882299108276066E-6, -4.7008654187390946E-6, -3.8724265256859586E-6,
      -3.2166445094620004E-6, -2.7304759035629134E-6, -2.3976739559322936E-6,
      -2.1924703946561381E-6, -2.083680499176413E-6, -2.0393800014851306E-6, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2650718801466323,
      0.58050741752112356, 0.845870876735914, 1.0307476451910988,
      1.1406251640786276, 1.1946555526458815, 1.2133525447025828,
      1.2131793539946374, 1.2052140053486344, 1.1958071146036526,
      1.1879185343050855, 1.1767917297989658, 0.8838387441031087,
      0.33511434377718163, -0.2330788898538764, -0.68150821173990694,
      -0.97933891282745578, -1.1476592785389421, -1.2243607727842492,
      -1.2457825991658213, -1.2394434524479196, -1.2230098099049291,
      -1.2060845194204346, -1.1516991541912385, -1.0438893677661507,
      -0.89361549145646, -0.72129630155850133, -0.54687575068642125,
      -0.38539593447497189, -0.24598200745847712, -0.13258210846213989,
      -0.045357563108382938, 0.0179030553152683, 0.060630755029117571,
      0.086659705809129017, 0.099716440365609407, 0.10314307701458608,
      0.099779559257655914, 0.091943084355267268, 0.081460424563014652,
      0.069724645459555729, 0.05775960948793233, 0.046283538290254805,
      0.035767606740377995, 0.026488059365479488, 0.018571540240434784,
      0.012033827806000848, 0.0068123482123254551, 0.0027929082776576076,
      -0.00016886926144963826, -0.0022308841866726858, -0.0035521080835337919,
      -0.0042837552381740137, -0.0045634797822301943, -0.0045120206907087483,
      -0.0042317544541721937, -0.0038066698646957355, -0.0033033442914509356,
      -0.0027725701090147762, -0.0022513485446717723, -0.0017650323794351384,
      -0.0013294562311203652, -0.00095294227511735687, -0.00063810983361552761,
      -0.0003834495519931122, -0.00018464757051395747, -3.5663118570783069E-5,
      7.0424662993624982E-5, 0.00014077709476514348, 0.00018233581182977546,
      0.00020151299994965353, 0.00020400353187496108, 0.0001946924938416465,
      0.00017763440672341639, 0.00015608369583464802, 0.00013255934099700795,
      0.00010892992285304522, 8.6508322883082066E-5, 6.6148041607547937E-5,
      4.833542600572976E-5, 3.3274031518438115E-5, 2.095889901201518E-5,
      1.1239732225275945E-5, 3.87285610574908E-6, -1.4375340680931367E-6,
      -5.0079121318653839E-6, -7.1525901824951866E-6, -8.1682483851793665E-6,
      -8.32453221098302E-6, -7.8594468728260056E-6, -6.9783151326676966E-6,
      -5.8550896960602223E-6, -4.6348168322998108E-6, -3.4360540274213377E-6,
      -2.3521170610658804E-6, -1.4503003990975309E-6, -7.68885431864861E-7,
      -3.1309899829380697E-7, -5.3500353745795536E-8, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0602875205865292, 1.2617421494979646,
      1.0614538368591622, 0.739507073820739, 0.43951007555011573,
      0.2161215542690148, 0.074787968226805251, -0.0006927628317813301,
      -0.031861394584012229, -0.03762756297992681, -0.031554321194268681,
      -0.044507218024478934, -1.1718119427834286, -2.1948976013037083,
      -2.272772934524232, -1.7937172875441223, -1.1913228043501953,
      -0.673281462845945, -0.30680597698122841, -0.085687305526288546,
      0.025356586871607131, 0.065734570171962145, 0.067701161937977788,
      0.21754146091678511, 0.43123914570035049, 0.60109550523876321,
      0.6892767595918351, 0.69768220348832, 0.64591926484579743,
      0.557655708065979, 0.45359959598534894, 0.34889818141502782,
      0.25304247369460492, 0.17091079885539709, 0.1041158031200458,
      0.052226938225921532, 0.01370654659590672, -0.013454071027720709,
      -0.03134589960955457, -0.041930639169010422, -0.0469431164138357,
      -0.047860143886493581, -0.045904284790710084, -0.042063726199507254,
      -0.037118189499594036, -0.031666076500178809, -0.026150849737735755,
      -0.020885918374701572, -0.016077759738671388, -0.011847110156428984,
      -0.0082480597008921883, -0.0052848955874444242, -0.0029265886185608856,
      -0.0011188981762247216, 0.00020583636608578393, 0.0011210649461462156,
      0.0017003383579058339, 0.0020133022929792006, 0.0021230967297446376,
      0.0020848862573720176, 0.0019452646609465351, 0.0017423045932590934,
      0.0015060558240120334, 0.0012593297660073173, 0.0010186411264896617,
      0.000795207925916619, 0.00059593780777269763, 0.00042435112625763217,
      0.00028140972708607406, 0.00016623486825852787, 7.6708752479512211E-5,
      9.9621277012301519E-6, -3.7244152133258293E-5, -6.8232348472920435E-5,
      -8.62028435550736E-5, -9.4097419350560208E-5, -9.4517672575850892E-5,
      -8.9686399879852573E-5, -8.1441125102136558E-5, -7.12504624072727E-5,
      -6.0245577949166589E-5, -4.9260530025691738E-5, -3.887666714695694E-5,
      -2.9467504478107462E-5, -2.1241560695368863E-5, -1.4281512255088988E-5,
      -8.5787122025192123E-6, -4.0626328107367188E-6, -6.2513530321460893E-7,
      1.8603413526280545E-6, 3.5245269606332362E-6, 4.4929017464298989E-6,
      4.8810914550416458E-6, 4.7950512195138917E-6, 4.3357478654218293E-6,
      3.6072666478733977E-6, 2.7256598689306797E-6, 1.8231457342842165E-6,
      1.0383945781920456E-6, 4.7908406832091486E-7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;

    helicopter_DW.From1_PWORK.TimePtr = (void *) pTimeValues0;
    helicopter_DW.From1_PWORK.DataPtr = (void *) pDataValues0;
    helicopter_DW.From1_IWORK.PrevIndex = 0;
  }

  /* Start for ToFile: '<Root>/To File' */
  {
    FILE *fp = (NULL);
    char fileName[509] = "states.mat";
    if ((fp = fopen(fileName, "wb")) == (NULL)) {
      rtmSetErrorStatus(helicopter_M, "Error creating .mat file states.mat");
      return;
    }

    if (rt_WriteMat4FileHeader(fp,8,0,"ans")) {
      rtmSetErrorStatus(helicopter_M,
                        "Error writing mat file header to file states.mat");
      return;
    }

    helicopter_DW.ToFile_IWORK.Count = 0;
    helicopter_DW.ToFile_IWORK.Decimation = -1;
    helicopter_DW.ToFile_PWORK.FilePtr = fp;
  }

  /* InitializeConditions for TransferFcn: '<S4>/Travel: Transfer Fcn' */
  helicopter_X.TravelTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Pitch: Transfer Fcn' */
  helicopter_X.PitchTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for TransferFcn: '<S4>/Elevation: Transfer Fcn' */
  helicopter_X.ElevationTransferFcn_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S7>/Integrator' */
  helicopter_X.Integrator_CSTATE = helicopter_P.Integrator_IC;

  /* InitializeConditions for Derivative: '<S4>/Derivative' */
  helicopter_DW.TimeStampA = (rtInf);
  helicopter_DW.TimeStampB = (rtInf);
}

/* Model terminate function */
void helicopter_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<Root>/HIL Initialize' */

  /* S-Function Block: helicopter/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_pwm_outputs = 0;
    hil_task_stop_all(helicopter_DW.HILInitialize_Card);
    hil_monitor_stop_all(helicopter_DW.HILInitialize_Card);
    is_switching = false;
    if ((helicopter_P.HILInitialize_set_analog_out_ex && !is_switching) ||
        (helicopter_P.HILInitialize_set_analog_outp_c && is_switching)) {
      {
        int_T i1;
        real_T *dw_AOVoltages = &helicopter_DW.HILInitialize_AOVoltages[0];
        for (i1=0; i1 < 8; i1++) {
          dw_AOVoltages[i1] = helicopter_P.HILInitialize_final_analog_outp;
        }
      }

      num_final_analog_outputs = 8U;
    }

    if ((helicopter_P.HILInitialize_set_pwm_output_ap && !is_switching) ||
        (helicopter_P.HILInitialize_set_pwm_outputs_p && is_switching)) {
      {
        int_T i1;
        real_T *dw_POValues = &helicopter_DW.HILInitialize_POValues[0];
        for (i1=0; i1 < 8; i1++) {
          dw_POValues[i1] = helicopter_P.HILInitialize_final_pwm_outputs;
        }
      }

      num_final_pwm_outputs = 8U;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_pwm_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(helicopter_DW.HILInitialize_Card
                         , helicopter_P.HILInitialize_analog_output_cha,
                         num_final_analog_outputs
                         , helicopter_P.HILInitialize_pwm_channels,
                         num_final_pwm_outputs
                         , NULL, 0
                         , NULL, 0
                         , &helicopter_DW.HILInitialize_AOVoltages[0]
                         , &helicopter_DW.HILInitialize_POValues[0]
                         , (t_boolean *) NULL
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog(helicopter_DW.HILInitialize_Card,
            helicopter_P.HILInitialize_analog_output_cha,
            num_final_analog_outputs, &helicopter_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_pwm_outputs > 0) {
          local_result = hil_write_pwm(helicopter_DW.HILInitialize_Card,
            helicopter_P.HILInitialize_pwm_channels, num_final_pwm_outputs,
            &helicopter_DW.HILInitialize_POValues[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(helicopter_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(helicopter_DW.HILInitialize_Card);
    hil_monitor_delete_all(helicopter_DW.HILInitialize_Card);
    hil_close(helicopter_DW.HILInitialize_Card);
    helicopter_DW.HILInitialize_Card = NULL;
  }

  /* Terminate for ToFile: '<Root>/To File' */
  {
    FILE *fp = (FILE *) helicopter_DW.ToFile_PWORK.FilePtr;
    if (fp != (NULL)) {
      char fileName[509] = "states.mat";
      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_M, "Error closing MAT-file states.mat");
        return;
      }

      if ((fp = fopen(fileName, "r+b")) == (NULL)) {
        rtmSetErrorStatus(helicopter_M, "Error reopening MAT-file states.mat");
        return;
      }

      if (rt_WriteMat4FileHeader(fp, 8, helicopter_DW.ToFile_IWORK.Count, "ans"))
      {
        rtmSetErrorStatus(helicopter_M,
                          "Error writing header for ans to MAT-file states.mat");
      }

      if (fclose(fp) == EOF) {
        rtmSetErrorStatus(helicopter_M, "Error closing MAT-file states.mat");
        return;
      }

      helicopter_DW.ToFile_PWORK.FilePtr = (NULL);
    }
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/

/* Solver interface called by GRT_Main */
#ifndef USE_GENERATED_SOLVER

void rt_ODECreateIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEDestroyIntegrationData(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

void rt_ODEUpdateContinuousStates(RTWSolverInfo *si)
{
  UNUSED_PARAMETER(si);
  return;
}                                      /* do nothing */

#endif

void MdlOutputs(int_T tid)
{
  helicopter_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  helicopter_update();
  UNUSED_PARAMETER(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  helicopter_initialize();
}

void MdlTerminate(void)
{
  helicopter_terminate();
}

/* Registration function */
RT_MODEL_helicopter_T *helicopter(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  helicopter_P.Integrator_UpperSat = rtInf;
  helicopter_P.Integrator_LowerSat = rtMinusInf;

  /* initialize real-time model */
  (void) memset((void *)helicopter_M, 0,
                sizeof(RT_MODEL_helicopter_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&helicopter_M->solverInfo,
                          &helicopter_M->Timing.simTimeStep);
    rtsiSetTPtr(&helicopter_M->solverInfo, &rtmGetTPtr(helicopter_M));
    rtsiSetStepSizePtr(&helicopter_M->solverInfo,
                       &helicopter_M->Timing.stepSize0);
    rtsiSetdXPtr(&helicopter_M->solverInfo, &helicopter_M->ModelData.derivs);
    rtsiSetContStatesPtr(&helicopter_M->solverInfo, (real_T **)
                         &helicopter_M->ModelData.contStates);
    rtsiSetNumContStatesPtr(&helicopter_M->solverInfo,
      &helicopter_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&helicopter_M->solverInfo,
      &helicopter_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&helicopter_M->solverInfo,
      &helicopter_M->ModelData.periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&helicopter_M->solverInfo,
      &helicopter_M->ModelData.periodicContStateRanges);
    rtsiSetErrorStatusPtr(&helicopter_M->solverInfo, (&rtmGetErrorStatus
      (helicopter_M)));
    rtsiSetRTModelPtr(&helicopter_M->solverInfo, helicopter_M);
  }

  rtsiSetSimTimeStep(&helicopter_M->solverInfo, MAJOR_TIME_STEP);
  helicopter_M->ModelData.intgData.f[0] = helicopter_M->ModelData.odeF[0];
  helicopter_M->ModelData.contStates = ((real_T *) &helicopter_X);
  rtsiSetSolverData(&helicopter_M->solverInfo, (void *)
                    &helicopter_M->ModelData.intgData);
  rtsiSetSolverName(&helicopter_M->solverInfo,"ode1");

  /* Initialize timing info */
  {
    int_T *mdlTsMap = helicopter_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    helicopter_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    helicopter_M->Timing.sampleTimes = (&helicopter_M->Timing.sampleTimesArray[0]);
    helicopter_M->Timing.offsetTimes = (&helicopter_M->Timing.offsetTimesArray[0]);

    /* task periods */
    helicopter_M->Timing.sampleTimes[0] = (0.0);
    helicopter_M->Timing.sampleTimes[1] = (0.002);

    /* task offsets */
    helicopter_M->Timing.offsetTimes[0] = (0.0);
    helicopter_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(helicopter_M, &helicopter_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = helicopter_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    helicopter_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(helicopter_M, -1);
  helicopter_M->Timing.stepSize0 = 0.002;
  helicopter_M->Timing.stepSize1 = 0.002;

  /* External mode info */
  helicopter_M->Sizes.checksums[0] = (3621657909U);
  helicopter_M->Sizes.checksums[1] = (3777463859U);
  helicopter_M->Sizes.checksums[2] = (4259787433U);
  helicopter_M->Sizes.checksums[3] = (3206432471U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    helicopter_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(helicopter_M->extModeInfo,
      &helicopter_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(helicopter_M->extModeInfo, helicopter_M->Sizes.checksums);
    rteiSetTPtr(helicopter_M->extModeInfo, rtmGetTPtr(helicopter_M));
  }

  helicopter_M->solverInfoPtr = (&helicopter_M->solverInfo);
  helicopter_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&helicopter_M->solverInfo, 0.002);
  rtsiSetSolverMode(&helicopter_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  helicopter_M->ModelData.blockIO = ((void *) &helicopter_B);

  {
    int32_T i;
    for (i = 0; i < 6; i++) {
      helicopter_B.Gain1[i] = 0.0;
    }

    helicopter_B.From = 0.0;
    helicopter_B.TravelCounttorad = 0.0;
    helicopter_B.Gain = 0.0;
    helicopter_B.Gain_d = 0.0;
    helicopter_B.PitchCounttorad = 0.0;
    helicopter_B.Gain_i = 0.0;
    helicopter_B.Gain_b = 0.0;
    helicopter_B.ElevationCounttorad = 0.0;
    helicopter_B.Gain_e = 0.0;
    helicopter_B.Sum = 0.0;
    helicopter_B.Gain_dg = 0.0;
    helicopter_B.Sum_g = 0.0;
    helicopter_B.From1[0] = 0.0;
    helicopter_B.From1[1] = 0.0;
    helicopter_B.From1[2] = 0.0;
    helicopter_B.From1[3] = 0.0;
    helicopter_B.Feedbackgain = 0.0;
    helicopter_B.Subtract1 = 0.0;
    helicopter_B.Sum1 = 0.0;
    helicopter_B.Sum2 = 0.0;
    helicopter_B.K_ei = 0.0;
    helicopter_B.Gain_l = 0.0;
    helicopter_B.BackmotorSaturation = 0.0;
    helicopter_B.FrontmotorSaturation = 0.0;
  }

  /* parameters */
  helicopter_M->ModelData.defaultParam = ((real_T *)&helicopter_P);

  /* states (continuous) */
  {
    real_T *x = (real_T *) &helicopter_X;
    helicopter_M->ModelData.contStates = (x);
    (void) memset((void *)&helicopter_X, 0,
                  sizeof(X_helicopter_T));
  }

  /* states (dwork) */
  helicopter_M->ModelData.dwork = ((void *) &helicopter_DW);
  (void) memset((void *)&helicopter_DW, 0,
                sizeof(DW_helicopter_T));

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_AIMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_AIMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_AOMinimums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_AOMaximums[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_AOVoltages[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_FilterFrequency[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_POSortedFreqs[i] = 0.0;
    }
  }

  {
    int32_T i;
    for (i = 0; i < 8; i++) {
      helicopter_DW.HILInitialize_POValues[i] = 0.0;
    }
  }

  helicopter_DW.TimeStampA = 0.0;
  helicopter_DW.LastUAtTimeA = 0.0;
  helicopter_DW.TimeStampB = 0.0;
  helicopter_DW.LastUAtTimeB = 0.0;
  helicopter_DW.HILWriteAnalog_Buffer[0] = 0.0;
  helicopter_DW.HILWriteAnalog_Buffer[1] = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    helicopter_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 16;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.B = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.P = &rtPTransTable;
  }

  /* Initialize Sizes */
  helicopter_M->Sizes.numContStates = (4);/* Number of continuous states */
  helicopter_M->Sizes.numPeriodicContStates = (0);/* Number of periodic continuous states */
  helicopter_M->Sizes.numY = (0);      /* Number of model outputs */
  helicopter_M->Sizes.numU = (0);      /* Number of model inputs */
  helicopter_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  helicopter_M->Sizes.numSampTimes = (2);/* Number of sample times */
  helicopter_M->Sizes.numBlocks = (67);/* Number of blocks */
  helicopter_M->Sizes.numBlockIO = (22);/* Number of block outputs */
  helicopter_M->Sizes.numBlockPrms = (148);/* Sum of parameter "widths" */
  return helicopter_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
