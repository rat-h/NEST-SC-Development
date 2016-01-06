/*
 *  stbrst_gc_conv.h
 *
 *  This file is part of a model to study development cortex's inputs 
 *  alignment in SC.
 * 
 *  Retinal module provides Starburst Amacrine cell model and
 *  Cholinergic connections modeled by convolution of voltage neighbor cells
 *
 *  Copyright (C) Ruben Tikidji-Hamburyan (rath@gwu.edu)
 *
 * * * * * * * * * * * *
 * 
 *  This module is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This module is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this module.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef STBRST_GC_CONV_H
#define STBRST_GC_CONV_H

#include "config.h"

#ifdef HAVE_GSL

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"

#include "universal_data_logger.h"
#include "recordables_map.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf_exp.h>

#include "convolution.h"
namespace nest
{
using std::vector;

/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 */
extern "C" int stbrst_gc_conv_dynamics( double, const double*, double*, void* );

/* BeginDocumentation
Name: stbrst_gc_conv - Starburst Amacrine cell model with convolution over 
 voltage neighbor Amacrine cells.

Description:

 stbrst_gc_conv is an implementation of a Starburst Amacrine cell with 
 extra of ganglion spiking neuron using the Hennig's model.
 The implementation use secondary connection event (same as gap-junction)
 to calculate voltage convolution for Cholinergic connections.

 This module is based on original code used for Henning et.al. JNeurosci. 2009
 kindly provided by Prof. Henning.

Parameters:
 In Hemming's model all parameters and variables is in SI. 
 We use scaled units for minimization of an rounding error:
 #  Cm   *  Vm   /  t    =  I   =   g   *   Vm
 # 1e-6F * 1e-3V / 1e-3s = 1e-6A = 1e-3S * 1e-3V
 #   uF  *  mV   /  ms   =  uA   =  mS   *  mV

 The following parameters can be set in the status dictionary.
 ========= STATE VARIABLES =========
 Vm         double - Membrane potential in mV 
 conCa      double - Cytosol Ca++ Concentration (unit ?)
 fAHP       double - Fast After-Hyperpolarization
 sAHP       double - Slow component of After-Hyperpolarization
 spont      double - Spontaneous activity, dynamical variable in uA
 synconv    double - Synaptic activity, dynamical variable in mS
 =========   AUX VARIABLE  =========
 VmThC		double - Thresholded Membrane positional for Convolution connections
                   - VmThC is Vm - VThC if Vm > VThC and the zero otherwise
 =========    PARAMETERS   =========
 VThC		double - Threshold Voltage Cholinergic activation
 
 El         double - Resting membrane potential in mV.
 gl         double - Leak conductance in mS.
 Cm         double - Capacity of the membrane in uF.
 
 
 gnoise     double - Conductance of spontaneous event in mS
 tnoise     double - Time constant of spontaneous integration in ms
 rnoise     double - Rate of spontaneous events in ms-1
 
 gsyn       double - Maximal conductance of synapses event in mS
 tsyn       double - Time constant of synaptic integration in ms
 
 eCa        double - Ca++ reversal potential in mV
 eK         double - K reversal potential in mV
 
 aFAHP      double - Fast AHP activation constant in (units ?).
 bFAHP      double - Fast AHP inactivation rate in ms-1.
 aSAHP      double - Slow AHP activation constant in (units ?).
 bSAHP      double - Slow AHP inactivation rate in ms-1.
 gAHP       double - Maximal conductance of cumulated AHP in mS
 
 vHCa       double - Half activation of Calcium steady-state in mV 
 vSCa       double - Slop of Calcium steady-state activation in mV-1
 gCa        double - Maximal conductance of Ca++ channel in mS
 gainCAcon  double - Gain of Cytosol Ca++ concentration in (unit ?)/uA
 tCAcon     double - Decay time constant of Ca++ concentration in ms
 
 totAHPinit double - * Initial condition for cumulative AHP (fAHP)
                     * In original model fAHP initiates from uniform random
                     *  distribution in range (0, 1). Set this parameter
                     *  individually for each neuron.
 
References: 
 Hennig, MH, Adams, C., Willshaw, D. and Sernagor, E. (2009).
 Early-stage retinal waves arise close to a critical state between
 local and global functional network connectivity.
 Journal of Neuroscience, 29: 1077-1086.

Sends: ConvolvEvent

Receives: SpikeEvent, ConvolvEvent, CurrentEvent, DataLoggingRequest

Author: Ruben Tikidji-Hamburyan
*/

namespace names{
	//States
	extern const Name Vm;
	extern const Name conCa;
	extern const Name fAHP;
	extern const Name sAHP;
	extern const Name spont;
	extern const Name synconv;
	//Aux Variable
	extern const Name VmThC;
	//Parameters
	extern const Name El;
	extern const Name gl;
	extern const Name Cm;
	extern const Name gnoise;
	extern const Name tnoise;
	extern const Name rnoise;
	extern const Name gsyn;
	extern const Name tsyn;
	extern const Name eCa;
	extern const Name eK;
	extern const Name aFAHP;
	extern const Name bFAHP;
	extern const Name aSAHP;
	extern const Name bSAHP;
	extern const Name gAHP;
	extern const Name vHCa;
	extern const Name vSCa;
	extern const Name gCa;
	extern const Name gainCAcon;
	extern const Name tCAcon;
	extern const Name totAHPinit;
}

class stbrst_gc_conv : public Archiving_Node
{

public:
  typedef Node base;

  stbrst_gc_conv();
  stbrst_gc_conv( const stbrst_gc_conv& );
  ~stbrst_gc_conv();

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and Hiding
   */
  using Node::handle;
  using Node::handles_test_event;

  port send_test_event( Node& target, rport receptor_type, synindex, bool );


  void handle( SpikeEvent& );
  void handle( CurrentEvent& );
  void handle( DataLoggingRequest& );
  void handle( ConvolvEvent& );

  port handles_test_event( SpikeEvent&, rport );
  port handles_test_event( CurrentEvent&, rport );
  port handles_test_event( DataLoggingRequest&, rport );
  port handles_test_event( ConvolvEvent&, rport );

  void
  sends_secondary_event( ConvolvEvent& )
  {
  }

  /**
   * Return membrane potential at time t.
potentials_.connect_logging_device();
   * This function is not thread-safe and should not be used in threaded
   * contexts to access the current membrane potential values.
   * @param Time the current network time
   *
   */
  double_t get_potential( Time const& ) const;

  /**
   * Define current membrane potential.
   * This function is thread-safe and should be used in threaded
   * contexts to change the current membrane potential value.
   * @param Time     the current network time
   * @param double_t new value of the mebrane potential
   *
   */
  void set_potential( Time const&, double_t );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

private:
  void init_state_( const Node& proto );
  void init_buffers_();
  void calibrate();
  bool update_( Time const&, const long_t, const long_t, const bool );
  void update( Time const&, const long_t, const long_t );
  bool prelim_update( Time const&, const long_t, const long_t );

  // END Boilerplate function declarations ----------------------------

  // Friends --------------------------------------------------------

  // make dynamics function quasi-member
  friend int stbrst_gc_conv_dynamics( double, const double*, double*, void* );

  // The next two classes need to be friend to access the State_ class/member
  friend class RecordablesMap< stbrst_gc_conv >;
  friend class UniversalDataLogger< stbrst_gc_conv >;

private:
  // ----------------------------------------------------------------

  //! Independent parameters
  struct Parameters_
  {
	//See documentation above :)
    double_t El_ ;
    double_t gl_ ;
    double_t Cm_ ;
    double_t gnoise_ ;
    double_t tnoise_ ;
    double_t rnoise_ ;
    double_t gsyn_ ;
    double_t tsyn_ ;
    double_t eCa_ ;
    double_t eK_ ;
    double_t aFAHP_ ;
    double_t bFAHP_ ;
    double_t aSAHP_ ;
    double_t bSAHP_ ;
    double_t gAHP_ ;
    double_t vHCa_ ;
    double_t vSCa_ ;
    double_t gCa_ ;
	double_t gainCAcon_ ;
	double_t tCAcon_ ;
	double_t totAHPinit_;

    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
    void set( const DictionaryDatum& ); //!< Set values from dicitonary
  };

public:
  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   * @note Copy constructor and assignment operator required because
   *       of C-style array.
   */
  struct State_
  {

    /**
     * Enumeration identifying elements in state array State_::y_.
     * The state vector must be passed to GSL as a C array. This enum
     * identifies the elements of the vector. It must be public to be
     * accessible from the iteration function.
     */
    enum StateVecElems
    {
      VM = 0,// 0
      CON_CA, // 1
      FAHP,   // 2
      SAHP,   // 3
	  SPONT,  // 4
	  SYNCONV,// 5
      STATE_VEC_SIZE
    };


    double_t y_[ STATE_VEC_SIZE ]; //!< neuron state, must be C-array for GSL solver
    int_t r_;                      //!< number of refractory steps remaining

    State_( const Parameters_& ); //!< Default initialization
    State_( const State_& );
    State_& operator=( const State_& );

    void get( DictionaryDatum& ) const;
    void set( const DictionaryDatum& );
  };

  // ----------------------------------------------------------------

private:
  /**
   * Buffers of the model.
   */
  struct Buffers_
  {
    Buffers_( stbrst_gc_conv& );                  //!<Sets buffer pointers to 0
    Buffers_( const Buffers_&, stbrst_gc_conv& ); //!<Sets buffer pointers to 0

    //! Logger for all analog data
    UniversalDataLogger< stbrst_gc_conv > logger_;

    /** buffers and sums up incoming spikes/currents */
    RingBuffer spike_exc_;
    RingBuffer spike_inh_;
    RingBuffer currents_;

    /** GSL ODE stuff */
    gsl_odeiv_step* s_;    //!< stepping function
    gsl_odeiv_control* c_; //!< adaptive stepsize control function
    gsl_odeiv_evolve* e_;  //!< evolution function
    gsl_odeiv_system sys_; //!< struct describing system

    // IntergrationStep_ should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double_t step_;          //!< step size in ms
    double IntegrationStep_; //!< current integration time step, updated by GSL

    // remembers current lag for piecewise interpolation
    long_t lag_;
    // remembers y_values from last prelim_update
    std::vector< double_t > last_y_values;
    // summarized gap weight
    double_t sumj_g_ij_;
    // summarized coefficients of the interpolation polynomial
    std::vector< double_t > interpolation_coefficients;

    /**
     * Input current injected by CurrentEvent.
     * This variable is used to transport the current applied into the
     * _dynamics function computing the derivative of the state vector.
     * It must be a part of Buffers_, since it is initialized once before
     * the first simulation, but not modified before later Simulate calls.
     */
    double_t I_stim_;
  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
    /** initial value to normalise excitatory synaptic current */
    double_t PSCurrInit_E_;

    /** initial value to normalise inhibitory synaptic current */
    double_t PSCurrInit_I_;

    int_t RefractoryCounts_;
  };

  // Access functions for UniversalDataLogger -------------------------------

  //! Read out state vector elements, used by UniversalDataLogger
  template < State_::StateVecElems elem >
  double_t
  get_y_elem_() const
  {
    return S_.y_[ elem ];
  }

  // ----------------------------------------------------------------

  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;

  //! Mapping of recordables names to access functions
  static RecordablesMap< stbrst_gc_conv > recordablesMap_;
};

inline void
stbrst_gc_conv::update( Time const& origin, const long_t from, const long_t to )
{
  update_( origin, from, to, false );
}

inline bool
stbrst_gc_conv::prelim_update( Time const& origin, const long_t from, const long_t to )
{
  bool done = false;
  State_ old_state = S_; // save state before prelim update
  done = update_( origin, from, to, true );
  S_ = old_state; // restore old state

  return done;
}

inline port
stbrst_gc_conv::send_test_event( Node& target, rport receptor_type, synindex, bool )
{
  SpikeEvent se;
  se.set_sender( *this );
  return target.handles_test_event( se, receptor_type );
}


inline port
stbrst_gc_conv::handles_test_event( SpikeEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
    throw UnknownReceptorType( receptor_type, get_name() );
  return 0;
}

inline port
stbrst_gc_conv::handles_test_event( CurrentEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
    throw UnknownReceptorType( receptor_type, get_name() );
  return 0;
}

inline port
stbrst_gc_conv::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
{
  if ( receptor_type != 0 )
    throw UnknownReceptorType( receptor_type, get_name() );
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline port
stbrst_gc_conv::handles_test_event( ConvolvEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
    throw UnknownReceptorType( receptor_type, get_name() );
  return 0;
}

inline void
stbrst_gc_conv::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d );
  Archiving_Node::get_status( d );

  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

inline void
stbrst_gc_conv::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  ptmp.set( d );         // throws if BadProperty
  State_ stmp = S_;      // temporary copy in case of errors
  stmp.set( d );         // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  Archiving_Node::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
}

} // namespace

#endif // HAVE_GSL
#endif // HH_PSC_ALPHA_GAP_H
