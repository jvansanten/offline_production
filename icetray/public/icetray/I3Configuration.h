/**
 *  $Id: I3Configuration.h 168750 2016-06-01 18:26:57Z cweaver $
 *  
 *  Copyright (C) 2007
 *  Troy D. Straszheim  <troy@icecube.umd.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *  
 *  This file is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>
 *  
 */
#ifndef ICETRAY_I3CONFIGURATION_H_INCLUDED
#define ICETRAY_I3CONFIGURATION_H_INCLUDED

#ifndef __CINT__
#include <string>
#include <map>
#include <iosfwd>
#include <boost/scoped_ptr.hpp>
#include <icetray/IcetrayFwd.h>
#include <icetray/I3DefaultName.h>
#include <icetray/serialization.h>
#include <icetray/is_shared_ptr.h>
#include <boost/python/extract.hpp>

/**
 * @brief This class holds the configuration.
 * 
 * @version $Id: I3Configuration.h 168750 2016-06-01 18:26:57Z cweaver $
 */

class I3ConfigurationImpl;

class I3Configuration
{
  boost::scoped_ptr<I3ConfigurationImpl> impl_;

public:

  virtual ~I3Configuration();

  I3Configuration();
  I3Configuration(const I3Configuration&);
  I3Configuration& operator = (const I3Configuration&);

  void Set(const std::string& key, const boost::python::object& value);

  bool Has(const std::string& name) const;

  // routine that does the talking to the impl
  void 
  Add(const std::string& name, const std::string& description, 
      const boost::python::object& default_value);

  // this one for w/o default value parameters
  void 
  Add(const std::string& name, const std::string& description); 

  template <typename T>
  void
  Add(const std::string& name, const std::string& description, 
      const T& default_value)
  {
    // forward to the boost::python::object version
    boost::python::object obj(default_value);
    Add(name, description, obj);
  };

  boost::python::object
  Get(const std::string& name) const;

  std::string
  GetDescription(const std::string& name) const;

  /**
   * This allows NoneType on the python side 
   * to translate to a Null pointer on the C++ side 
   */
  template <typename T>
    typename boost::enable_if<is_shared_ptr<T>, T>::type
  Get(const std::string& name) const
  {
    boost::python::object obj(Get(name));
    if(obj.ptr() == Py_None) return T();
    return boost::python::extract<T>(obj);
  }

  template <typename T>
    typename boost::disable_if<is_shared_ptr<T>, T>::type
  Get(const std::string& name) const
  {
    boost::python::object obj(Get(name));
    return boost::python::extract<T>(obj);
  }

  std::string ClassName() const;
  void ClassName(const std::string& s);

  std::string InstanceName() const;
  void InstanceName(const std::string& s);

  bool is_ok() const;

  ///
  /// prints in xml format for consumption by xsltproc.  You can't just use
  /// the icecube::serialization xml since this isn't well formed
  /// due to <px class_id="9" class_name="typeholder<bool>" ...
  ///
  std::string inspect() const;

  std::vector<std::string> keys() const;

private:

  friend class icecube::serialization::access;

  template <typename Archive>
  void
  serialize(Archive&, unsigned);

  friend std::ostream& operator<<(std::ostream&, const I3Configuration&);

};

std::ostream& operator<<(std::ostream&, const I3Configuration&);

I3_DEFAULT_NAME(I3Configuration);
I3_POINTER_TYPEDEFS(I3Configuration);

#endif // ndef __CINT__

#endif // I3CONFIGURATION_H
