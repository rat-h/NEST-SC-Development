/*
 *  convolution.h
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


/* BeginDocumentation
Name: gap_junction - Synapse type for gap-junction connections.

Description:
 gap_junction is a connector to create gap junctions between pairs
 of neurons. Please note that gap junctions are two-way connections:
 In order to create an accurate gap-junction connection between two
 neurons i and j two connections are required:

 i j conn_spec gap_junction   Connect
 j i conn_spec gap_junction   Connect

 The value of the parameter "delay" is ignored for connections of
 type gap_junction.

Transmits: GapJEvent

References:

 Hahne, J., Helias, M., Kunkel, S., Igarashi, J.,
 Bolten, M., Frommer, A. and Diesmann, M.,
 A unified framework for spiking and gap-junction interactions
 in distributed neuronal network simulations,
 Front. Neuroinform. 9:22. (2015),
 doi: 10.3389/fninf.2015.00022

 Mancilla, J. G., Lewis, T. J., Pinto, D. J.,
 Rinzel, J., and Connors, B. W.,
 Synchronization of electrically coupled pairs
 of inhibitory interneurons in neocortex,
 J. Neurosci. 27, 2058-2073 (2007),
 doi: 10.1523/JNEUROSCI.2715-06.2007

Author: Jan Hahne, Moritz Helias, Susanne Kunkel
SeeAlso: synapsedict, hh_psc_alpha_gap
*/


#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include "connection.h"

namespace nest
{
class DSConvolvEvent;
/**
 * ConvolvEvent is used for send voltage over 
 * Retinal to make a convolution.
 */
class ConvolvEvent : public Event
{
  double_t v_;

public:
  void operator()();
  ConvolvEvent* clone() const;
  void event_hook( DSConvolvEvent& );

  void set_voltage( double_t );
  double_t get_voltage() const;
};

inline ConvolvEvent*
ConvolvEvent::clone() const
{
  return new ConvolvEvent( *this );
}

inline void
ConvolvEvent::set_voltage( double_t v )
{
  v_ = v;
}

inline double_t
ConvolvEvent::get_voltage() const
{
  return v_;
}

/**
 * 'Callback request event'
 */
class DSConvolvEvent : public ConvolvEvent
{
public:
  void operator()();
};



//template < typename targetidentifierT >
//class Convolv : public Connection< targetidentifierT >
template < typename targetidentifierT >
class Convolv : public nest::Connection< targetidentifierT >
{

  double_t weight_;

public:
  // this line determines which common properties to use
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;
  //typedef ConvolvEvent EventType;


  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  Convolv()
    : ConnectionBase()
    , weight_( 1.0 )
  {
  }


  // Explicitly declare all methods inherited from the dependent base ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not automatically
  // found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;
  class ConnTestDummyNode : public nest::ConnTestDummyNodeBase
  {
  public:
    using nest::ConnTestDummyNodeBase::handles_test_event;
    nest::port
    handles_test_event( nest::ConvolvEvent&, nest::rport )
    {
      return nest::invalid_port_;
    }

    nest::port
    handles_test_event( nest::DSConvolvEvent&, nest::rport )
    {
      return nest::invalid_port_;
    }
  };


  void
  check_connection( Node& s, Node& t, rport receptor_type, double_t, const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;
    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );
  }

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param p The port under which this connection is stored in the Connector.
   * \param t_lastspike Time point of last spike emitted
   */
  void
  send( Event& e, thread t, double_t, const CommonSynapseProperties& )
  {
    e.set_weight( weight_ );
    //e.set_delay( ConnectionBase::get_delay_steps() );
    //e.set_receiver( *ConnectionBase::get_target( t ) );
    //e.set_rport( ConnectionBase::get_rport() );
    e.set_delay( get_delay_steps() );
    e.set_receiver( *get_target( t ) );
    e.set_rport( get_rport() );
    e();
  }

  void get_status( DictionaryDatum& d ) const;

  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  void
  set_weight( double_t w )
  {
    weight_ = w;
  }
};

template < typename targetidentifierT >
void
Convolv< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double_t >( d, names::weight, weight_ );
  def< long_t >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
Convolv< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double_t >( d, names::weight, weight_ );
}

} // namespace

#endif /* #ifndef GAP_JUNCTION_H */
