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
#include "convolution.h"

nest::RecordablesMap< nest::stbrst_gc_conv > nest::stbrst_gc_conv::recordablesMap_;
namespace nest
{
//Names exported to the NEST

void ConvolvEvent::operator()()
{
  ((stbrst_gc_conv*) receiver_)->handle( *this );
}


void DSConvolvEvent::operator()()
{
 ((stbrst_gc_conv*)sender_)->event_hook( *this );
}

// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< stbrst_gc_conv >::create()
{
  // use standard names whereever you can for consistency!
  insert_( "Vm", &stbrst_gc_conv::get_y_elem_< stbrst_gc_conv::State_::VM > );
  insert_( "conCa", &stbrst_gc_conv::get_y_elem_< stbrst_gc_conv::State_::CON_CA > );
  insert_( "fAHP", &stbrst_gc_conv::get_y_elem_< stbrst_gc_conv::State_::FAHP > );
  insert_( "sAHP", &stbrst_gc_conv::get_y_elem_< stbrst_gc_conv::State_::SAHP > );
  insert_( "spont", &stbrst_gc_conv::get_y_elem_< stbrst_gc_conv::State_::SPONT > );
  insert_( "synconv", &stbrst_gc_conv::get_y_elem_< stbrst_gc_conv::State_::SYNCONV > );
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
  f[S::VM     ] = 
      (IsAHP + ICa + spont - synconv * V + P.gl_ * (P.El_ - V) ) / P.Cm_;
  f[S::CON_CA ] = P.gainCAcon_* ICa - conCa / P.tCAcon_;
  f[S::SPONT  ] = - spont / P.tnoise_;
  f[S::FAHP   ] = (P.aFAHP_ * conCa + sAHP) * (1. - fAHP) - fAHP * P.bFAHP_;
  f[S::SAHP   ] = P.aSAHP_ * conCa4/(2e-10 + conCa4) * (1. - sAHP) - sAHP * P.bSAHP_;
  f[S::SYNCONV] = - synconv / P.tsyn_;


  return GSL_SUCCESS;
}

} //nest namespace


/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::stbrst_gc_conv::Parameters_::Parameters_()
	: VThC_ ( -60. ) //in mV
    , El_( -65. ) // in mV
    , gl_( 1e3/180000000 ) // in mS
    , Cm_( 1.6E-4 ) // in uF
    , gnoise_( 2E-7 ) // in mS
    , tnoise_( 80. ) // in ms
    , rnoise_( 1.4 ) // in 1/ms
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
    , vSCa_( 10. ) // in mV-1
    , gCa_( 1e-5 ) // in mS // -!
	, gainCAcon_ ( 1e2 ) // in ? / uA
	, tCAcon_ ( 50. ) // in ms
	, totAHPinit_ (0.5) // in ?
	, seed_ (0) // in ?
{
}

nest::stbrst_gc_conv::State_::State_( const Parameters_& p)
  : r_( 0 )
{
  const double_t vinit = (y_[ 0 ] = p.El_);

  for ( size_t i = 1; i < STATE_VEC_SIZE; ++i )
    y_[ i ] = 0;

  y_[ FAHP ] = p.totAHPinit_;
  //DB>>
  //printf("\nRTH_DB: Set y_[FAHP]=%g\n",p.totAHPinit_);
  //<<DB
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
	def< double >( d, "VThC", VThC_ );
	def< double >( d, "El", El_ );
	def< double >( d, "gl", gl_ );
	def< double >( d, "Cm", Cm_ );
	def< double >( d, "gnoise", gnoise_ );
	def< double >( d, "tnoise", tnoise_ );
	def< double >( d, "rnoise", rnoise_ );
	def< double >( d, "gsyn", gsyn_ );
	def< double >( d, "tsyn", tsyn_ );
	def< double >( d, "eCa", eCa_ );
	def< double >( d, "eK", eK_ );
	def< double >( d, "aFAHP", aFAHP_ );
	def< double >( d, "bFAHP", bFAHP_ );
	def< double >( d, "aSAHP", aSAHP_ );
	def< double >( d, "bSAHP", bSAHP_ );
	def< double >( d, "gAHP", gAHP_ );
	def< double >( d, "vHCa", vHCa_ );
	def< double >( d, "vSCa", vSCa_ );
	def< double >( d, "gCa", gCa_ );
	def< double >( d, "gainCAcon", gainCAcon_ );
	def< double >( d, "tCAcon", tCAcon_ );
	def< double >( d, "totAHPinit", totAHPinit_ );
	def< long_t >( d, "seed", seed_ );
}

void
nest::stbrst_gc_conv::Parameters_::set( const DictionaryDatum& d )
{
	 updateValue< double >( d, "VThC", VThC_ );
	 updateValue< double >( d, "El", El_ );
	 updateValue< double >( d, "gl", gl_ );
	 updateValue< double >( d, "Cm", Cm_ );
	 updateValue< double >( d, "gnoise", gnoise_ );
	 updateValue< double >( d, "tnoise", tnoise_ );
	 updateValue< double >( d, "rnoise", rnoise_ );
	 updateValue< double >( d, "gsyn", gsyn_ );
	 updateValue< double >( d, "tsyn", tsyn_ );
	 updateValue< double >( d, "eCa", eCa_ );
	 updateValue< double >( d, "eK", eK_ );
	 updateValue< double >( d, "aFAHP", aFAHP_ );
	 updateValue< double >( d, "bFAHP", bFAHP_ );
	 updateValue< double >( d, "aSAHP", aSAHP_ );
	 updateValue< double >( d, "bSAHP", bSAHP_ );
	 updateValue< double >( d, "gAHP", gAHP_ );
	 updateValue< double >( d, "vHCa", vHCa_ );
	 updateValue< double >( d, "vSCa", vSCa_ );
	 updateValue< double >( d, "gCa", gCa_ );
	 updateValue< double >( d, "gainCAcon", gainCAcon_ );
	 updateValue< double >( d, "tCAcon", tCAcon_ );
	 updateValue< double >( d, "totAHPinit", totAHPinit_ );
	 updateValue< long_t >( d, "seed", seed_ );
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
	def< double >( d, "Vm", y_[ VM ] );
	def< double >( d, "conCa", y_[ CON_CA ] );
	def< double >( d, "fAHP", y_[ FAHP ] );
	def< double >( d, "sAHP", y_[ SAHP ] );
	def< double >( d, "spont", y_[ SPONT ] );
	def< double >( d, "synconv", y_[ SYNCONV ] );
}

void
nest::stbrst_gc_conv::State_::set( const DictionaryDatum& d )
{
	updateValue< double >( d, "Vm", y_[ VM ] );
	updateValue< double >( d, "conCa", y_[ CON_CA ] );
	updateValue< double >( d, "fAHP", y_[ FAHP ] );
	updateValue< double >( d, "sAHP", y_[ SAHP ] );
	updateValue< double >( d, "spont", y_[ SPONT ] );
	updateValue< double >( d, "synconv", y_[ SYNCONV ] );
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
  , gslrnd( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest::stbrst_gc_conv::Buffers_::Buffers_( const Buffers_&, stbrst_gc_conv& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
  , gslrnd( 0 )
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
//  Node::set_needs_prelim_update( true );
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
  if ( B_.gslrnd )
    gsl_rng_free( B_.gslrnd );
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
  B_.syn_convol_.clear(); // includes resize

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
  /**
   * Int random number generator
   */
  B_.gslrnd = gsl_rng_alloc( gsl_rng_default );
  gsl_rng_set(B_.gslrnd, ( P_.seed_ ) ? P_.seed_ : get_thread() );
}

void
nest::stbrst_gc_conv::calibrate()
{
  B_.logger_.init(); // ensures initialization in case mm connected after Simulate
  S_.y_[ State_::FAHP ] = P_.totAHPinit_;
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
nest::stbrst_gc_conv::update( Time const& origin,
  const long_t from,
  const long_t to )
{

  assert( to >= 0 && ( delay ) from < Scheduler::get_min_delay() );
  assert( from < to );


  for ( long_t lag = from; lag < to; ++lag )
  {

    double tt = 0.0; // it's all relative!
    const double V = S_.y_[ State_::VM ];
    const double_t ca_act   = (0.5 * (1. + std::tanh((V - P_.vHCa_) / (2. * P_.vSCa_))));

	// First, we accumulate all voltage event to avoid delay = min_delay + 1step 
    S_.y_[ State_::SYNCONV ] += B_.syn_convol_.get_value( lag ) * P_.gsyn_;
	

	// Original model has 
	if (gsl_rng_uniform(B_.gslrnd) < (1. - std::exp(-P_.rnoise_ * ca_act * (1. - ca_act) * B_.step_)))
        S_.y_[ State_::SPONT ] += P_.gnoise_ * (P_.eCa_ - V);

    // adaptive step integration
    while ( tt < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply( B_.e_,
        B_.c_,
        B_.s_,
        &B_.sys_,             // system of ODE
        &tt,                  // from t...
        B_.step_,             // ...to t=t+h
        &B_.IntegrationStep_, // integration window (written on!)
        S_.y_ );              // neuron state


      if ( status != GSL_SUCCESS )
        throw GSLSolverFailure( get_name(), status );
    }

    
	// GC spike
	if ( S_.y_[ State_::VM ] >= P_.VThC_ && V < P_.VThC_ ) {
		set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

		SpikeEvent se;
		network()->send( *this, se, lag );
	}
    
	if ( S_.y_[ State_::VM ] > P_.VThC_ ){
		//S_.r_ = V_.RefractoryCounts_;

		//set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

		ConvolvEvent se;
		se.set_voltage(S_.y_[ State_::VM ] - P_.VThC_);
		network()->send( *this, se, lag );
	}

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }

  return;
}

void
nest::stbrst_gc_conv::handle( SpikeEvent& e )
{
}


void
nest::stbrst_gc_conv::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

void
nest::stbrst_gc_conv::handle( ConvolvEvent& e )
{
  assert( e.get_delay() > 0 );
  const double_t v = e.get_voltage();
  const double_t w = e.get_weight();
  B_.syn_convol_.add_value( e.get_rel_delivery_steps( network()->get_slice_origin() ), w * v );
}

#endif // HAVE_GSL
