/*
 *  stbrst_gc_conv.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "stbrst_gc_conv.h"

#ifdef HAVE_GSL

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include <limits>
#include "universal_data_logger_impl.h"

#include <iomanip>
#include <iostream>
#include <cstdio>
#include <cmath> // in case we need isnan() // fabs
#include <nest_names.h>

nest::RecordablesMap< nest::stbrst_gc_conv > nest::stbrst_gc_conv::recordablesMap_;
namespace nest
{
//Names exported to the NEST
namespace names{
	const Name Vm( "Vm" );
	const Name conCa( "conCa" );
	const Name fAHP( "fAHP" );
	const Name sAHP( "sAHP" );
	const Name spont( "spont" );
	const Name synconv( "synconv" );
	const Name VmThC( "VmThC" );
	const Name El( "El" );
	const Name gl( "gl" );
	const Name Cm( "Cm" );
	const Name gnoise( "gnoise" );
	const Name tnoise( "tnoise" );
	const Name rnoise( "rnoise" );
	const Name gsyn( "gsyn" );
	const Name tsyn( "tsyn" );
	const Name eCa( "eCa" );
	const Name eK( "eK" );
	const Name aFAHP( "aFAHP" );
	const Name bFAHP( "bFAHP" );
	const Name aSAHP( "aSAHP" );
	const Name bSAHP( "bSAHP" );
	const Name gAHP( "gAHP" );
	const Name vHCa( "vHCa" );
	const Name vSCa( "vSCa" );
	const Name gCa( "gCa" );
	const Name gainCAcon( "gainCAcon" );
	const Name tCAcon( "tCAcon" );
	const Name totAHPinit( "totAHPinit" );
}

// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< stbrst_gc_conv >::create()
{
  // use standard names whereever you can for consistency!
  insert_( names::Vm, &stbrst_gc_conv::get_y_elem_< stbrst_gc_conv::State_::VM > );
  insert_( names::conCa, &stbrst_gc_conv::get_y_elem_< stbrst_gc_conv::State_::CON_CA > );
  insert_( names::fAHP, &stbrst_gc_conv::get_y_elem_< stbrst_gc_conv::State_::FAHP > );
  insert_( names::sAHP, &stbrst_gc_conv::get_y_elem_< stbrst_gc_conv::State_::SAHP > );
  insert_( names::spont, &stbrst_gc_conv::get_y_elem_< stbrst_gc_conv::State_::SPONT > );
  insert_( names::synconv, &stbrst_gc_conv::get_y_elem_< stbrst_gc_conv::State_::SYNCONV > );
}

extern "C" int
stbrst_gc_conv_dynamics( double time, const double y[], double f[], void* pnode )
{
  // a shorthand
  typedef nest::stbrst_gc_conv::State_ S;

  // get access to node so we can almost work as in a member function
  assert( pnode );
  const nest::stbrst_gc_conv& node = *( reinterpret_cast< nest::stbrst_gc_conv* >( pnode ) );
  
  //parameters shorthand
  const nest::stbrst_gc_conv::Parameters_ P = node.P_;
  
  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[].

  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...

  // shorthand for state variables
  const double_t& V       = y[ S::VM ];
  const double_t& conCa   = y[ S::CON_CA ];
  const double_t& fAHP    = y[ S::FAHP ];
  const double_t& sAHP    = y[ S::SAHP ];
  const double_t& spont   = y[ S::SPONT ];
  const double_t& synconv = y[ S::SYNCONV ];
  
  const double_t ca_act   = (0.5 * (1. + std::tanh((V - P.vHCa_) / (2. * P.vSCa_))));
  const double_t ICa      = ca_act * P.gCa_ * (P.eCa_ - V);
  const double_t dc       = fAHP*fAHP*fAHP*fAHP;
  const double_t IsAHP    = P.gAHP_ * dc * (P.eK_ - V);
  const double_t conCa4   = conCa * conCa * conCa * conCa;
  f[   S::VM  ] =
      (IsAHP + ICa + spont - synconv * V + P.gl_ * (P.El_ - V) ) / P.Cm_;
  f[ S::SPONT ] = - spont / P.tnoise_;
  f[ S::FAHP  ] = P.aFAHP_ * conCa + sAHP * (1. - fAHP) - fAHP * P.bFAHP_;
  f[ S::SAHP  ] = P.aSAHP_ * conCa4/(2e-10 + conCa4) * (1. - sAHP) - sAHP * P.bSAHP_;



/*


  // set I_gap depending on interpolation order
  double_t gap = 0.0;

  const double_t t = time / node.B_.step_;

  switch ( Scheduler::get_prelim_interpolation_order() )
  {
  case 0:
    gap = -node.B_.sumj_g_ij_ * V + node.B_.interpolation_coefficients[ node.B_.lag_ ];
    break;

  case 1:
    gap = -node.B_.sumj_g_ij_ * V + node.B_.interpolation_coefficients[ node.B_.lag_ * 2 + 0 ]
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 2 + 1 ] * t;
    break;

  case 3:
    gap = -node.B_.sumj_g_ij_ * V + node.B_.interpolation_coefficients[ node.B_.lag_ * 4 + 0 ]
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 4 + 1 ] * t
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 4 + 2 ] * t * t
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 4 + 3 ] * t * t * t;
    break;

  default:
    throw BadProperty( "Interpolation order must be 0, 1, or 3." );
  }

  const double_t I_gap = gap;

  // V dot -- synaptic input are currents, inhib current is negative
  f[ S::VM ] =
    ( -( I_Na + I_K + I_L ) + node.B_.I_stim_ + node.P_.I_e + I_ex + I_in + I_gap ) / node.P_.C_m;

  // channel dynamics
  f[ S::HH_M ] = alpha_m * ( 1 - y[ S::HH_M ] ) - beta_m * y[ S::HH_M ]; // m-variable
  f[ S::HH_H ] = alpha_h * ( 1 - y[ S::HH_H ] ) - beta_h * y[ S::HH_H ]; // h-variable
  f[ S::HH_P ] = alpha_p * ( 1 - y[ S::HH_P ] ) - beta_p * y[ S::HH_P ]; // p-variable
  f[ S::HH_N ] = alpha_n * ( 1 - y[ S::HH_N ] ) - beta_n * y[ S::HH_N ]; // n-variable

  // synapses: alpha functions
  f[ S::DI_EXC ] = -dI_ex / node.P_.tau_synE;
  f[ S::I_EXC ] = dI_ex - ( I_ex / node.P_.tau_synE );
  f[ S::DI_INH ] = -dI_in / node.P_.tau_synI;
  f[ S::I_INH ] = dI_in - ( I_in / node.P_.tau_synI );
*/
  return GSL_SUCCESS;
}
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::stbrst_gc_conv::Parameters_::Parameters_()
    : El_( -65. ) // in mV
    , gl_( 1e3/180000000 ) // in mS
    , Cm_( 1.6E-4 ) // in uF
    , gnoise_( 2E-7 ) // in mS
    , tnoise_( 80. ) // in ms
    , gsyn_( 3.2e-4 ) // in mS
    , tsyn_( 200. ) // in ms
    , eCa_( 50. ) // in mV
    , eK_( -90. ) // in mV
    , aFAHP_( 2400. ) // in ?
    , bFAHP_( 200. ) // in mc
    , aSAHP_( 30. ) // in ?
    , bSAHP_( 20. ) // in mc
    , gAHP_( 3e-5 ) // in mS // -(minus) !!!!
    , vHCa_( -20. ) // in mV
    , vSCa_( 1. ) // in mV-1
    , gCa_( 1e-5 ) // in mS // -!
	, gainCAcon_ ( 1e2 ) // in ? / uA
	, tCAcon_ ( 50. ) // in ms
	, totAHPinit_ (0.5) // in ?
{
}

nest::stbrst_gc_conv::State_::State_( const Parameters_& p)
  : r_( 0 )
{
  const double_t vinit = (y_[ 0 ] = p.El_);

  for ( size_t i = 1; i < STATE_VEC_SIZE; ++i )
    y_[ i ] = 0;

  y_[ FAHP ] = p.totAHPinit_;
}

nest::stbrst_gc_conv::State_::State_( const State_& s )
  : r_( s.r_ )
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
    y_[ i ] = s.y_[ i ];
}

nest::stbrst_gc_conv::State_& nest::stbrst_gc_conv::State_::operator=( const State_& s )
{
  assert( this != &s ); // would be bad logical error in program

  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
    y_[ i ] = s.y_[ i ];
  r_ = s.r_;
  return *this;
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::stbrst_gc_conv::Parameters_::get( DictionaryDatum& d ) const
{
	def< double >( d, names::El, El_ );
	def< double >( d, names::gl, gl_ );
	def< double >( d, names::Cm, Cm_ );
	def< double >( d, names::gnoise, gnoise_ );
	def< double >( d, names::tnoise, tnoise_ );
	def< double >( d, names::gsyn, gsyn_ );
	def< double >( d, names::tsyn, tsyn_ );
	def< double >( d, names::eCa, eCa_ );
	def< double >( d, names::eK, eK_ );
	def< double >( d, names::aFAHP, aFAHP_ );
	def< double >( d, names::bFAHP, bFAHP_ );
	def< double >( d, names::aSAHP, aSAHP_ );
	def< double >( d, names::bSAHP, bSAHP_ );
	def< double >( d, names::gAHP, gAHP_ );
	def< double >( d, names::vHCa, vHCa_ );
	def< double >( d, names::vSCa, vSCa_ );
	def< double >( d, names::gCa, gCa_ );
	def< double >( d, names::gainCAcon, gainCAcon_ );
	def< double >( d, names::tCAcon, tCAcon_ );
	def< double >( d, names::totAHPinit, totAHPinit_ );
}

void
nest::stbrst_gc_conv::Parameters_::set( const DictionaryDatum& d )
{
	 updateValue< double >( d, names::El, El_ );
	 updateValue< double >( d, names::gl, gl_ );
	 updateValue< double >( d, names::Cm, Cm_ );
	 updateValue< double >( d, names::gnoise, gnoise_ );
	 updateValue< double >( d, names::tnoise, tnoise_ );
	 updateValue< double >( d, names::gsyn, gsyn_ );
	 updateValue< double >( d, names::tsyn, tsyn_ );
	 updateValue< double >( d, names::eCa, eCa_ );
	 updateValue< double >( d, names::eK, eK_ );
	 updateValue< double >( d, names::aFAHP, aFAHP_ );
	 updateValue< double >( d, names::bFAHP, bFAHP_ );
	 updateValue< double >( d, names::aSAHP, aSAHP_ );
	 updateValue< double >( d, names::bSAHP, bSAHP_ );
	 updateValue< double >( d, names::gAHP, gAHP_ );
	 updateValue< double >( d, names::vHCa, vHCa_ );
	 updateValue< double >( d, names::vSCa, vSCa_ );
	 updateValue< double >( d, names::gCa, gCa_ );
	 updateValue< double >( d, names::gainCAcon, gainCAcon_ );
	 updateValue< double >( d, names::tCAcon, tCAcon_ );
	 updateValue< double >( d, names::totAHPinit, totAHPinit_ );
	if ( Cm_ <= 0 )
		throw BadProperty( "Capacitance must be strictly positive." );
////// TODO >> Check parameters range (rth)
/*  if ( t_ref_ < 0 )
    throw BadProperty( "Refractory time cannot be negative." );

  if ( tau_synE <= 0 || tau_synI <= 0 )
    throw BadProperty( "All time constants must be strictly positive." );

  if ( g_Kv1 < 0 || g_Kv3 < 0 || g_Na < 0 || g_L < 0 )
    throw BadProperty( "All conductances must be non-negative." );
*/
}

void
nest::stbrst_gc_conv::State_::get( DictionaryDatum& d ) const
{
	def< double >( d, names::Vm, y_[ VM ] );
	def< double >( d, names::conCa, y_[ CON_CA ] );
	def< double >( d, names::fAHP, y_[ FAHP ] );
	def< double >( d, names::sAHP, y_[ SAHP ] );
	def< double >( d, names::spont, y_[ SPONT ] );
	def< double >( d, names::synconv, y_[ SYNCONV ] );
}

void
nest::stbrst_gc_conv::State_::set( const DictionaryDatum& d )
{
	updateValue< double >( d, names::Vm, y_[ VM ] );
	updateValue< double >( d, names::conCa, y_[ CON_CA ] );
	updateValue< double >( d, names::fAHP, y_[ FAHP ] );
	updateValue< double >( d, names::sAHP, y_[ SAHP ] );
	updateValue< double >( d, names::spont, y_[ SPONT ] );
	updateValue< double >( d, names::synconv, y_[ SYNCONV ] );
////// TODO >> Check parameters range (rth)
/*	if ( y_[ HH_M ] < 0 || y_[ HH_H ] < 0 || y_[ HH_N ] < 0 || y_[ HH_P ] < 0 )
		throw BadProperty( "All (in)activation variables must be non-negative." );
*/
}

nest::stbrst_gc_conv::Buffers_::Buffers_( stbrst_gc_conv& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest::stbrst_gc_conv::Buffers_::Buffers_( const Buffers_&, stbrst_gc_conv& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::stbrst_gc_conv::stbrst_gc_conv()
  : Archiving_Node()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  Node::set_needs_prelim_update( true );
  recordablesMap_.create();
}

nest::stbrst_gc_conv::stbrst_gc_conv( const stbrst_gc_conv& n )
  : Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

nest::stbrst_gc_conv::~stbrst_gc_conv()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ )
    gsl_odeiv_step_free( B_.s_ );
  if ( B_.c_ )
    gsl_odeiv_control_free( B_.c_ );
  if ( B_.e_ )
    gsl_odeiv_evolve_free( B_.e_ );
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest::stbrst_gc_conv::init_state_( const Node& proto )
{
  const stbrst_gc_conv& pr = downcast< stbrst_gc_conv >( proto );
  S_ = pr.S_;
}

void
nest::stbrst_gc_conv::init_buffers_()
{
  B_.spike_exc_.clear(); // includes resize
  B_.spike_inh_.clear(); // includes resize
  B_.currents_.clear();  // includes resize

  // allocate strucure for gap events here
  // function is called from Scheduler::prepare_nodes() before the
  // first call to update
  // so we already know which interpolation scheme to use according
  // to the properties of this neurons
  // determine size of structure depending on interpolation scheme
  // and unsigned int Scheduler::min_delay() (number of simulation time steps per min_delay step)

  // resize interpolation_coefficients depending on interpolation order
  const size_t quantity =
    Scheduler::get_min_delay() * ( Scheduler::get_prelim_interpolation_order() + 1 );

  B_.interpolation_coefficients.resize( quantity, 0.0 );

  B_.last_y_values.resize( Scheduler::get_min_delay(), 0.0 );

  B_.sumj_g_ij_ = 0.0;

  Archiving_Node::clear_history();

  B_.logger_.reset();

  B_.step_ = Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE );
  else
    gsl_odeiv_step_reset( B_.s_ );

  if ( B_.c_ == 0 )
    B_.c_ = gsl_odeiv_control_y_new( 1e-6, 0.0 );
  else
    gsl_odeiv_control_init( B_.c_, 1e-6, 0.0, 1.0, 0.0 );

  if ( B_.e_ == 0 )
    B_.e_ = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  else
    gsl_odeiv_evolve_reset( B_.e_ );

  B_.sys_.function = stbrst_gc_conv_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );

  B_.I_stim_ = 0.0;
}

void
nest::stbrst_gc_conv::calibrate()
{
/*  B_.logger_.init(); // ensures initialization in case mm connected after Simulate

  V_.PSCurrInit_E_ = 1.0 * numerics::e / P_.tau_synE;
  V_.PSCurrInit_I_ = 1.0 * numerics::e / P_.tau_synI;
  V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref_ ) ).get_steps();
  assert( V_.RefractoryCounts_ >= 0 ); // since t_ref_ >= 0, this can only fail in error
*/}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

bool
nest::stbrst_gc_conv::update_( Time const& origin,
  const long_t from,
  const long_t to,
  const bool prelim )
{

  assert( to >= 0 && ( delay ) from < Scheduler::get_min_delay() );
  assert( from < to );

  bool done = true;
  const size_t interpolation_order = Scheduler::get_prelim_interpolation_order();
  const double_t prelim_tol = Scheduler::get_prelim_tol();

  // allocate memory to store the new interpolation coefficients
  // to be sent by gap event
  const size_t quantity = Scheduler::get_min_delay() * ( interpolation_order + 1 );
  std::vector< double_t > new_coefficients( quantity, 0.0 );

  // parameters needed for piecewise interpolation
  double_t y_i = 0.0, y_ip1 = 0.0, hf_i = 0.0, hf_ip1 = 0.0;
  double_t f_temp[ State_::STATE_VEC_SIZE ];

  for ( long_t lag = from; lag < to; ++lag )
  {

    // B_.lag is needed by stbrst_gc_conv_dynamics to
    // determine the current section
    B_.lag_ = lag;

    if ( prelim ){ //Preliminary update
      y_i = S_.y_[ State_::VM ];
      if ( interpolation_order == 3 )
      {
        stbrst_gc_conv_dynamics( 0, S_.y_, f_temp, reinterpret_cast< void* >( this ) );
        hf_i = B_.step_ * f_temp[ State_::VM ];
      }
    }

    double_t t = 0.0;
    const double_t U_old = S_.y_[ State_::VM ];

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( t < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply( B_.e_,
        B_.c_,
        B_.s_,
        &B_.sys_,             // system of ODE
        &t,                   // from t
        B_.step_,             // to t <= step
        &B_.IntegrationStep_, // integration step size
        S_.y_ );              // neuronal state

      if ( status != GSL_SUCCESS )
        throw GSLSolverFailure( get_name(), status );
    }

    if ( !prelim ){ // main update
      //S_.y_[ State_::DI_EXC ] += B_.spike_exc_.get_value( lag ) * V_.PSCurrInit_E_;
      //S_.y_[ State_::DI_INH ] += B_.spike_inh_.get_value( lag ) * V_.PSCurrInit_I_;
      // sending spikes: crossing 0 mV, pseudo-refractoriness and local maximum...
      // refractory?
      if ( S_.r_ > 0 )
        --S_.r_;
      else
        // (    threshold    &&     maximum       )
        if ( S_.y_[ State_::VM ] >= 0 && U_old > S_.y_[ State_::VM ] )
      {
        S_.r_ = V_.RefractoryCounts_;

        set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

        SpikeEvent se;
        network()->send( *this, se, lag );
      }

      // log state data
      B_.logger_.record_data( origin.get_steps() + lag );

      // set new input current
      B_.I_stim_ = B_.currents_.get_value( lag );
    }
    else // if(prelim)
    {
      //S_.y_[ State_::DI_EXC ] += B_.spike_exc_.get_value_prelim( lag ) * V_.PSCurrInit_E_;
      //S_.y_[ State_::DI_INH ] += B_.spike_inh_.get_value_prelim( lag ) * V_.PSCurrInit_I_;
      // check deviation from last iteration
      done = ( fabs( S_.y_[ State_::VM ] - B_.last_y_values[ lag ] ) <= prelim_tol ) && done;
      B_.last_y_values[ lag ] = S_.y_[ State_::VM ];

      // update different interpolations

      // constant term is the same for each interpolation order
      new_coefficients[ lag * ( interpolation_order + 1 ) + 0 ] = y_i;

      switch ( interpolation_order )
      {
      case 0:
        break;

      case 1:
        y_ip1 = S_.y_[ State_::VM ];

        new_coefficients[ lag * ( interpolation_order + 1 ) + 1 ] = y_ip1 - y_i;
        break;

      case 3:
        y_ip1 = S_.y_[ State_::VM ];
        stbrst_gc_conv_dynamics( B_.step_, S_.y_, f_temp, reinterpret_cast< void* >( this ) );
        hf_ip1 = B_.step_ * f_temp[ State_::VM ];

        new_coefficients[ lag * ( interpolation_order + 1 ) + 1 ] = hf_i;
        new_coefficients[ lag * ( interpolation_order + 1 ) + 2 ] =
          -3 * y_i + 3 * y_ip1 - 2 * hf_i - hf_ip1;
        new_coefficients[ lag * ( interpolation_order + 1 ) + 3 ] =
          2 * y_i - 2 * y_ip1 + hf_i + hf_ip1;
        break;

      default:
        throw BadProperty( "Interpolation order must be 0, 1, or 3." );
      }
    }


  } // end for-loop

  // if !prelim perform constant extrapolation and reset last_y_values
  if ( !prelim ) //main update
  {
    for ( long_t temp = from; temp < to; ++temp )
      new_coefficients[ temp * ( interpolation_order + 1 ) + 0 ] = S_.y_[ State_::VM ];

    B_.last_y_values.clear();
    B_.last_y_values.resize( Scheduler::get_min_delay(), 0.0 );
  }

  // Send gap-event
  ConvolvEvent ge;
  ge.set_coeffarray( new_coefficients );
  network()->send_secondary( *this, ge );

  // Reset variables
  B_.sumj_g_ij_ = 0.0;
  B_.interpolation_coefficients.clear();
  B_.interpolation_coefficients.resize( quantity, 0.0 );

  return done;
}

void
nest::stbrst_gc_conv::handle( SpikeEvent& e )
{
  assert( e.get_delay() > 0 );

  if ( e.get_weight() > 0.0 )
    B_.spike_exc_.add_value( e.get_rel_delivery_steps( network()->get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  else
    B_.spike_inh_.add_value( e.get_rel_delivery_steps( network()->get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() ); // current input, keep negative weight
}

void
nest::stbrst_gc_conv::handle( CurrentEvent& e )
{
  assert( e.get_delay() > 0 );

  const double_t c = e.get_current();
  const double_t w = e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value( e.get_rel_delivery_steps( network()->get_slice_origin() ), w * c );
}

void
nest::stbrst_gc_conv::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

void
nest::stbrst_gc_conv::handle( ConvolvEvent& e )
{

  B_.sumj_g_ij_ += e.get_weight();

  size_t i = 0;
  std::vector< uint_t >::iterator it = e.begin();
  // The call to get_coeffvalue( it ) in this loop also advances the iterator it
  while ( it != e.end() )
  {
    B_.interpolation_coefficients[ i++ ] += e.get_weight() * e.get_coeffvalue( it );
  }
}

#endif // HAVE_GSL
