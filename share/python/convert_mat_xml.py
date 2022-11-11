#!/usr/bin/env python
from __future__ import division, print_function

import argparse
from lxml import etree
from lxml import objectify
import os.path
import sys

def _clean_ns(tag):
  if not isinstance(tag, str):
    return None
  if '}' in tag:
    tag = tag.split('}')[1]
  return tag

def _get_ns(tag):
  if not isinstance(tag, str):
    return None
  if '}' in tag:
    tag = tag.split('}')[0][1:]
  return tag

class mat_converter:
  
  nsmap_old = {None : "http://www.cfs++.org"}
  nsmap_new = {None : "http://www.cfs++.org/material"}

  def __init__(self):
    self.E = objectify.ElementMaker(annotate=False,
                                    namespace=self.nsmap_new[None],
                                    nsmap=self.nsmap_new)
  
  def convert_mat_xml(self, infile, outfile):
    if not os.path.exists(infile):
      print("Input file does not exist", file=sys.stderr, flush=True)
      return 1
    
    parser = etree.XMLParser(ns_clean=True, recover=True)
    tree = etree.parse(infile, parser)
    
    inroot = tree.getroot()
    if (_clean_ns(inroot.tag) != "cfsMaterialDataBase" or
        _get_ns(inroot.tag) != self.nsmap_old[None]):
      print("Input file is not an old openCFS material database", file=sys.stderr, flush=True)
      return 1
    
    outroot = self.E.cfsMaterialDataBase()
    
    for mat in inroot:
      if isinstance(mat, etree._Comment) > 0:
        outroot.append(etree.Comment(mat.text))
      else:
        outroot.append(self._convert_material(mat))
    
    with open(outfile, "wb") as out:
      out.write(etree.tostring(outroot, encoding="utf-8", pretty_print=True, xml_declaration=True))
    
    return 0
  
  def _ensure_elem(self, parent, elem):
    result = parent.find(elem)
    
    if result is None:
      result = etree.Element(elem, nsmap=self.nsmap_new)
      parent.append(result)
    
    return result
  
  def _create_scalar(self, name, real_value, imag_value=None):
    result = etree.Element(name, nsmap=self.nsmap_new)
    result.append(self.E.real(real_value))
    if imag_value is not None:
      result.append(self.E.imag(imag_value))
       
    return result
  
  def _create_linscalar(self, real_value, imag_value=None):
    return self._create_scalar("linear", real_value, imag_value)
  
  def _convert_linscalar(self, scalar, parent):
    result = etree.Element(_clean_ns(scalar.tag), nsmap=self.nsmap_new)
    result.append(self._create_linscalar(scalar.text))
    parent.append(result)
    
  def _create_tensor(self, name, real, imag=None, dim1=None, dim2=None):
    result = etree.Element(name, nsmap=self.nsmap_new)
    result.append(self.E.real(real))
    if imag is not None:
      result.append(self.E.imag(imag))
    if dim1 is not None:
      result.set("dim1", dim1)
    if dim2 is None:
      dim2 = dim1
    if dim2 is not None:
      result.set("dim2", dim2)
    return result
  
  def _create_lintensor(self, real, imag=None, dim1=None, dim2=None):
    return self._create_tensor("linear", real, imag, dim1, dim2)
  
  def _create_orthotensor(self, real1, real2, real3, imag1=None, imag2=None, imag3=None):
    return self.E.orthotropic(self._create_scalar("value_1", real1, imag1),
                              self._create_scalar("value_2", real2, imag2),
                              self._create_scalar("value_3", real3, imag3))
    
  def _convert_nonlintensor(self, prop, parent):
    lin = prop.find("./linear", namespaces=self.nsmap_old)
    if lin is not None:
      self._convert_tensorprop(lin, parent)
    
    nl = prop.find("./nonlinear", namespaces=self.nsmap_old)
    if nl is not None:
      self._convert_nonlin(nl, parent)
  
  def _convert_tensorprop(self, prop, parent):
    iso_real = prop.findtext("./isotropic/real", namespaces=self.nsmap_old)
    iso_imag = prop.findtext("./isotropic/imag", namespaces=self.nsmap_old)
    if iso_real is None:
      iso_real = prop.findtext("./isotropic", namespaces=self.nsmap_old)
    if iso_real is not None:
      parent.append(self.E.linear(self._create_scalar("isotropic", iso_real, iso_imag)))
    
    tensor = prop.find("./tensor", namespaces=self.nsmap_old)
    if tensor is not None:
      parent.append(self.E.linear(self._convert_tensor("tensor", tensor)))

  def _convert_tensor(self, name, tensor):
    real = tensor.findtext("./real", namespaces=self.nsmap_old)
    imag = tensor.findtext("./imag", namespaces=self.nsmap_old)
    if real is None:
      real = tensor.text
    if real is not None:
      dim1 = tensor.get("dim1")
      dim2 = tensor.get("dim2")
      return self._create_tensor(name, real, imag, dim1, dim2)
    
  def _convert_nonlin(self, nl, parent):
    result = self.E.nonlinear()
    
    iso = nl.find("./isotropic", namespaces=self.nsmap_old)
    transiso = nl.find("./transversalIsotropic", namespaces=self.nsmap_old)
    ortho = nl.find("./orthotropic", namespaces=self.nsmap_old)
    tensor = nl.find("./tensor", namespaces=self.nsmap_old)

    if iso is not None:
      new_iso = self.E.isotropic()
      result.append(new_iso)
      self._convert_nonlindata(new_iso, iso)
    if transiso is not None:
      new_transiso = self.E.transversalIsotropic()
      result.append(new_transiso)
      self._convert_nonlindata(new_transiso, transiso)
    if ortho is not None:
      new_ortho = self.E.orthotropic()
      result.append(new_ortho)
      self._convert_nonlindata(new_ortho, ortho)
    if tensor is not None:
      new_tensor = self.E.tensor()
      result.append(new_tensor)
      self._convert_nonlindata(new_tensor, tensor)

    parent.append(result)
      
  def _convert_nonlindata(self, result, nl):
    dep = nl.findtext("./dependency", namespaces=self.nsmap_old)
    approxType = nl.findtext("./approxType", namespaces=self.nsmap_old)
    accu = nl.findtext("./measAccuracy", namespaces=self.nsmap_old)
    approxVal = nl.findtext("./maxApproxVal", namespaces=self.nsmap_old)
    dataname = nl.findtext("./dataName", namespaces=self.nsmap_old)
    entry = nl.findtext("./entry", namespaces=self.nsmap_old)
    nlattr = nl.get("nonlinear")
    
    result.append(self.E.dependency(dep if dep is not None else ""))
    result.append(self.E.approxType(approxType if approxType is not None else "No approximation"))
    result.append(self.E.measAccuracy(accu if accu is not None else "0"))
    result.append(self.E.maxApproxVal(approxVal if approxVal is not None else "0"))
    result.append(self.E.dataName(dataname if dataname is not None else ""))
    
    if entry is not None:
      result.append(self.E.entry(entry))
    if nlattr is not None:
      result.attrib["nonlinear"] = nlattr
    
  def _deepcopy(self, elem):
    result = etree.Element(_clean_ns(elem.tag), nsmap=self.nsmap_new)
  
    if elem.text is not None and len(elem.text) > 0:
      result.text = elem.text
    
    for sub in elem:
      if isinstance(sub, etree._Comment):
        result.append(etree.Comment(sub.text))
      else:
        result.append(self._deepcopy(sub))
    
    for (attrname, attrval) in elem.items():
      result.set(attrname, attrval)
    
    return result
  
  def _convert_material(self, mat):
    result = self.E.material(name=mat.get('name'))
    
    # loop over sub-elements
    for pde in mat:
      if isinstance(pde, etree._Comment):
        result.append(etree.Comment(pde.text))
        continue
      
      tagname = _clean_ns(pde.tag)
      if tagname == "mechanical":
        self._convert_mechanical(pde, result)
      elif tagname == "acoustic":
        self._convert_acoustic(pde, result)
      elif tagname == "electric":
        self._convert_electric(pde, result)
      elif tagname == "flow":
        self._convert_flow(pde, result)
      elif tagname == "testmat":
        self._convert_testmat(pde, result)
      elif tagname == "heatConduction":
        self._convert_heat(pde, result)
      elif tagname == "elecConduction":
        self._convert_elecconduction(pde, result)
      elif tagname == "magnetic":
        self._convert_magnetic(pde, result)
      elif tagname == "magnetoStrictive":
        self._convert_magstrict(pde, result)
      elif tagname == "piezo":
        self._convert_piezo(pde, result)
      elif tagname == "pyroelectric":
        self._convert_pyro(pde, result)
      elif tagname == "thermoelastic":
        self._convert_thermoelast(pde, result)
    
    return result
  
  def _convert_density(self, dens, parent):
    parent.append(self.E.density(self._create_linscalar(dens.text)))
  
  def _convert_mechanical(self, mech, mat):
    result = self.E.mechanical()
    
    for prop in mech:
      if isinstance(prop, etree._Comment):
        result.append(etree.Comment(prop.text))
        continue
      
      tagname = _clean_ns(prop.tag)
      if tagname == "density":
        self._convert_density(prop, result)
      elif tagname == "elasticity":
        self._convert_elast(prop, result)
      elif tagname == "viscoelasticity":
        print("Warning: viscoelasticity properties were ignored", file=sys.stderr)
      elif tagname == "elasticityCoefficient":
        elast = self._ensure_elem(result, "elasticity")
        self._convert_nonlin(prop, elast)
      #elif tagname == "irreversibleStrainCoefficient":
      #  self._convert_irrstraincoef(prop, result)
      elif tagname == "thermalExpanison":
        self._convert_thermexp(prop, result)
      elif tagname == "mechanicalDamping":
        self._convert_mechdamping(prop, result)
      elif tagname == "magnetoStrictionTensor_h_mech":
        real = prop.findtext("./real", namespaces=self.nsmap_old)
        imag = prop.findtext("./imag", namespaces=self.nsmap_old)
        if real is not None:
          dim1 = prop.get("dim1")
          dim2 = prop.get("dim2")
          result.append(self._create_tensor("magnetoStrictionTensor_h_mech", real, imag, dim1, dim2))
    
    mat.append(result)
        
  def _convert_elast(self, elast, mech):
    result = self._ensure_elem(mech, "elasticity")
    lin = self._ensure_elem(result, "linear")

    oldType = elast.find("./*")
    matType = _clean_ns(oldType.tag)
    typeElem = self._ensure_elem(lin, matType)
    
    if matType == "isotropic":
      e_real = oldType.findtext("./real/elasticityModulus", namespaces=self.nsmap_old)
      e_imag = oldType.findtext("./imag/elasticityModulus", namespaces=self.nsmap_old)
      nu_real = oldType.findtext("./real/poissonNumber", namespaces=self.nsmap_old)
      nu_imag = oldType.findtext("./imag/poissonNumber", namespaces=self.nsmap_old)
      k_real = oldType.findtext("./real/compressionModulus", namespaces=self.nsmap_old)
      k_imag = oldType.findtext("./imag/compressionModulus", namespaces=self.nsmap_old)
      g_real = oldType.findtext("./real/shearModulus", namespaces=self.nsmap_old)
      g_imag = oldType.findtext("./imag/shearModulus", namespaces=self.nsmap_old)
      l_real = oldType.findtext("./real/lameParameterLamda", namespaces=self.nsmap_old)
      l_imag = oldType.findtext("./imag/lameParameterLamda", namespaces=self.nsmap_old)
      mu_real = oldType.findtext("./real/lameParameterMu", namespaces=self.nsmap_old)
      mu_imag = oldType.findtext("./imag/lameParameterMu", namespaces=self.nsmap_old)
      
      if e_real is not None and nu_real is not None:
        typeElem.append(self._create_scalar("elasticityModulus", e_real, e_imag))
        typeElem.append(self._create_scalar("poissonNumber", nu_real, nu_imag))
      elif k_real is not None and g_real is not None:
        typeElem.append(self._create_scalar("shearModulus", g_real, g_imag))
        typeElem.append(self._create_scalar("compressionModulus", k_real, k_imag))
      elif l_real is not None and mui_real is not None:
        typeElem.append(self._create_scalar("lameParameterMu", mu_real, mu_imag))
        typeElem.append(self._create_scalar("lameParameterLambda", l_real, l_imag))
    
    elif matType == "transversalIsotropic":
      e_real = oldType.findtext("./real/elasticityModulus", namespaces=self.nsmap_old)
      e_imag = oldType.findtext("./imag/elasticityModulus", namespaces=self.nsmap_old)
      e3_real = oldType.findtext("./real/elasticityModulus_3", namespaces=self.nsmap_old)
      e3_imag = oldType.findtext("./imag/elasticityModulus_3", namespaces=self.nsmap_old)
      nu_real = oldType.findtext("./real/poissonNumber", namespaces=self.nsmap_old)
      nu_imag = oldType.findtext("./imag/poissonNumber", namespaces=self.nsmap_old)
      nu3_real = oldType.findtext("./real/poissonNumber_3", namespaces=self.nsmap_old)
      nu3_imag = oldType.findtext("./imag/poissonNumber_3", namespaces=self.nsmap_old)
      g_real = oldType.findtext("./real/shearModulus", namespaces=self.nsmap_old)
      g_imag = oldType.findtext("./imag/shearModulus", namespaces=self.nsmap_old)
      g3_real = oldType.findtext("./real/shearModulus_3", namespaces=self.nsmap_old)
      g3_imag = oldType.findtext("./imag/shearModulus_3", namespaces=self.nsmap_old)
      
      typeElem.append(self._create_scalar("elasticityModulus", e_real, e_imag))
      typeElem.append(self._create_scalar("elasticityModulus_3", e3_real, e3_imag))
      typeElem.append(self._create_scalar("poissonNumber", nu_real, nu_imag))
      typeElem.append(self._create_scalar("poissonNumber_3", nu3_real, nu3_imag))
      typeElem.append(self._create_scalar("shearModulus", g_real, g_imag))
      typeElem.append(self._create_scalar("shearModulus_3", g3_real, g3_imag))
    
    elif matType == "orthotropic":
      e1_real = oldType.findtext("./real/elasticityModulus_1", namespaces=self.nsmap_old)
      e1_imag = oldType.findtext("./imag/elasticityModulus_1", namespaces=self.nsmap_old)
      e2_real = oldType.findtext("./real/elasticityModulus_2", namespaces=self.nsmap_old)
      e2_imag = oldType.findtext("./imag/elasticityModulus_2", namespaces=self.nsmap_old)
      e3_real = oldType.findtext("./real/elasticityModulus_3", namespaces=self.nsmap_old)
      e3_imag = oldType.findtext("./imag/elasticityModulus_3", namespaces=self.nsmap_old)
      nu12_real = oldType.findtext("./real/poissonNumber_12", namespaces=self.nsmap_old)
      nu12_imag = oldType.findtext("./imag/poissonNumber_12", namespaces=self.nsmap_old)
      nu23_real = oldType.findtext("./real/poissonNumber_23", namespaces=self.nsmap_old)
      nu23_imag = oldType.findtext("./imag/poissonNumber_23", namespaces=self.nsmap_old)
      nu13_real = oldType.findtext("./real/poissonNumber_13", namespaces=self.nsmap_old)
      nu13_imag = oldType.findtext("./imag/poissonNumber_13", namespaces=self.nsmap_old)
      g12_real = oldType.findtext("./real/shearModulus_12", namespaces=self.nsmap_old)
      g12_imag = oldType.findtext("./imag/shearModulus_12", namespaces=self.nsmap_old)
      g23_real = oldType.findtext("./real/shearModulus_23", namespaces=self.nsmap_old)
      g23_imag = oldType.findtext("./imag/shearModulus_23", namespaces=self.nsmap_old)
      g13_real = oldType.findtext("./real/shearModulus_13", namespaces=self.nsmap_old)
      g13_imag = oldType.findtext("./imag/shearModulus_13", namespaces=self.nsmap_old)
      
      typeElem.append(self._create_scalar("elasticityModulus_1", e1_real, e1_imag))
      typeElem.append(self._create_scalar("elasticityModulus_2", e2_real, e2_imag))
      typeElem.append(self._create_scalar("elasticityModulus_3", e3_real, e3_imag))
      typeElem.append(self._create_scalar("poissonNumber_12", nu12_real, nu12_imag))
      typeElem.append(self._create_scalar("poissonNumber_23", nu23_real, nu23_imag))
      typeElem.append(self._create_scalar("poissonNumber_13", nu13_real, nu13_imag))
      typeElem.append(self._create_scalar("shearModulus_12", g12_real, g12_imag))
      typeElem.append(self._create_scalar("shearModulus_23", g23_real, g23_imag))
      typeElem.append(self._create_scalar("shearModulus_13", g13_real, g13_imag))
      
    elif matType == "tensor":
      real_tensor = oldType.findtext("./real", namespaces=self.nsmap_old)
      imag_tensor = oldType.findtext("./imag", namespaces=self.nsmap_old)
      dim1 = oldType.get("dim1")
      dim2 = oldType.get("dim2") if "dim2" in oldType.keys() else dim1
      
      typeElem.set("dim1", dim1)
      typeElem.set("dim2", dim2)
      typeElem.append(self.E.real(real_tensor))
      if imag_tensor is not None:
        typeElem.append(self.E.imag(imag_tensor))
  
  def _convert_irrstraincoef(self, irr, parent):
    coeffs = irr.findtext("./coeffs", namespaces=self.nsmap_old)
    parent.append(self.E.coeffs(coeffs))
  
  def _convert_thermexp(self, exp, parent):
    result = self.E.thermalExpansion()
    
    iso_real = exp.findtext("./isotropic/real", namespaces=self.nsmap_old)
    iso_imag = exp.findtext("./isotropic/imag", namespaces=self.nsmap_old)
    ortho_real = exp.findtext("./orthotropic/real", namespaces=self.nsmap_old)
    ortho_imag = exp.findtext("./orthotropic/imag", namespaces=self.nsmap_old)
    aniso_real = exp.findtext("./anisotropic/real", namespaces=self.nsmap_old)
    aniso_imag = exp.findtext("./anisotropic/imag", namespaces=self.nsmap_old)
    temp_real = exp.findtext("./refTemperature/real", namespaces=self.nsmap_old)
    temp_imag = exp.findtext("./refTemperature/imag", namespaces=self.nsmap_old)
    
    if iso_real is not None:
      result.append(self._create_scalar("isotropic", iso_real, iso_imag))
    elif ortho_real is not None:
      or_vec = ortho_real.split()
      oi_vec = ortho_imag.split() if ortho_imag is not None else [None, None, None] 
      result.append(self._create_orthotensor(or_vec[0], or_vec[1], or_vec[2],
                                             oi_vec[0], oi_vec[1], oi_vec[2]))
    elif aniso_real is not None:
      result.append(self._create_tensor("tensor", aniso_real, aniso_imag, "6", "1"))
    
    if temp_real is not None:
      result.append(self._create_scalar("refTemperature", temp_real, temp_imag))
      
    parent.append(result)

  def _convert_mechdamping(self, damp, parent):
    result = self.E.damping()
    
    # Rayleigh damping
    self._convert_rayleighdamp(damp, result)
    
    ## fractional damping
    #alg = damp.findtext("./fractional/alg", namespaces=self.nsmap_old)
    #mem = damp.findtext("./fractional/memory", namespaces=self.nsmap_old)
    #interpol = damp.findtext("./fractional/interpolation", namespaces=self.nsmap_old)
    #
    #if alg is not None and mem is not None and interpol is not None:
    #  result.append(self.E.fractional(self.E.alg(alg),
    #                                  self.E.memory(mem),
    #                                  self.E.interpolation(interpol)))
    
    parent.append(result)

  def _convert_rayleighdamp(self, damp, parent):
    alpha = damp.findtext("./rayleigh/alpha", namespaces=self.nsmap_old)
    beta = damp.findtext("./rayleigh/beta", namespaces=self.nsmap_old)
    tanDelta = damp.findtext("./rayleigh/lossTangensDelta", namespaces=self.nsmap_old)
    freq = damp.findtext("./rayleigh/measuredFreq", namespaces=self.nsmap_old)
  
    if alpha is not None and beta is not None:
      parent.append(self.E.rayleigh(self.E.alpha(alpha),
                                    self.E.beta(beta),
                                    self.E.measuredFreq(freq)))
    elif tanDelta is not None:
      parent.append(self.E.rayleigh(self.E.lossTangensDelta(tanDelta),
                                    self.E.measuredFreq(freq)))
    
  def _convert_acoustic(self, acou, parent):
    result = self.E.acoustic()
  
    for prop in acou:
      if isinstance(prop, etree._Comment):
        result.append(etree.Comment(prop.text))
        continue
      
      tagname = _clean_ns(prop.tag)
      if tagname == "density":
        self._convert_density(prop, result)
      elif tagname == "densityComplex":
        dens_real = prop.findtext('./real', namespaces=self.nsmap_old)
        dens_imag = prop.findtext('./imag', namespaces=self.nsmap_old)
        if dens_real is not None:
          result.append(self.E.densityComplex(self._create_linscalar(dens_real, dens_imag)))
      elif (tagname == "compressionModulus" or
          tagname == "compressionModulusComplex"):
        self._convert_acoucompmod(prop, result)
      elif tagname == "acousticDamping":
        self._convert_acoudamping(prop, result)
      elif tagname == "acousticNonlinear":
        bovera = prop.findtext("./acousticNonlinear/bOverA", namespaces=self.nsmap_old)
        if bovera is not None:
          result.append(self.E.nonlinear(self._create_scalar("bOverA", bovera)))
      elif tagname == "kinematicViscosity":
        self._convert_linscalar(prop, result)
      elif tagname == "adiabaticExponent":
        self._convert_linscalar(prop, result)
    
    parent.append(result)
  
  def _convert_acoucompmod(self, k, parent):
    real = k.findtext("./real", namespaces=self.nsmap_old)
    imag = k.findtext("./imag", namespaces=self.nsmap_old)
    
    if real is None:
      real = k.text
    
    if imag is None:
      result = self.E.compressionModulus()
    else:
      result = self.E.compressionModulusComplex()
    
    result.append(self._create_linscalar(real, imag))
    parent.append(result)
  
  def _convert_acoudamping(self, damp, parent):
    result = self.E.damping()
    
    # Rayleigh damping
    self._convert_rayleighdamp(damp, result)
    
    ## thermo-viscous damping
    #alpha0 = damp.findtext("./thermoViscous/alpha0", namespaces=self.nsmap_old)
    #if alpha0 is not None:
    #  result.append(self._create_scalar("thermoViscous", alpha0))
    
    ## Fractional damping
    #alpha0 = damp.findtext("./fractional/alpha0", namespaces=self.nsmap_old)
    #y = damp.findtext("./fractional/y", namespaces=self.nsmap_old)
    #if alpha0 is not None and y is not None:
    #  result.append(self.E.fractional(self._create_scalar("alpha0", alpha0),
    #                                  self._create_scalar("y", y)))
    
    parent.append(result)

  def _convert_electric(self, elec, parent):
    result = self.E.electric()
  
    for prop in elec:
      if isinstance(prop, etree._Comment):
        result.append(etree.Comment(prop.text))
        continue
      
      tagname = _clean_ns(prop.tag)
      if tagname == "permittivity":
        permit = self._ensure_elem(result, "permittivity")
        self._convert_tensorprop(prop, permit)
      elif tagname == "permittivityCoefficient":
        permit = self._ensure_elem(result, "permittivity")
        self._convert_nonlin(prop, permit)
      elif tagname == "hystModel":
        result.append(self._deepcopy(prop))
    
    parent.append(result)
  
  def _convert_elecconduction(self, elec, parent):
    result = self.E.elecConduction()
    
    for prop in elec:
      if isinstance(prop, etree._Comment):
        result.append(etree.Comment(prop.text))
        continue
      
      if _clean_ns(prop.tag) == "elecConductivity":
        self._convert_elecconductivity(prop, result)
    
    parent.append(result)
    
  def _convert_elecconductivity(self, cond, parent):
    result = self.E.electricConductivity()
    
    # support really old syntax
    if len(cond) == 0:
      result.append(self.E.linear(self._create_scalar("isotropic", cond.text)))
    
    self._convert_nonlintensor(cond, result)
#     lin = cond.find("./linear")
#     if lin is not None:
#       newlin = self.E.linear()
#       self._convert_tensorprop(lin, newlin)
#       result.append(newlin)
#     
#     nl = cond.find("./nonlinear")
#     if nl is not None:
#       nonlin = self.E.nonlinear()
#       self._convert_nonlin(nl, result)
#       result.append(nonlin)
    
    parent.append(result)
    
  def _convert_flow(self, flow, parent):
    result = self.E.flow()
    
    for prop in flow:
      if isinstance(prop, etree._Comment):
        result.append(etree.Comment(prop.text))
        continue
      
      tagname = _clean_ns(prop.tag)
      if tagname == "density":
        self._convert_density(prop, result)
      elif tagname == "dynamicViscosity":
        result.append(self.E.dynamicViscosity(self._create_linscalar(prop.text)))
      elif tagname == "kinematicViscosity":
        result.append(self.E.kinematicViscosity(self._create_linscalar(prop.text)))
      elif tagname == "bulkViscosity":
        result.append(self.E.bulkViscosity(self._create_linscalar(prop.text)))
      elif tagname == "adiabaticExponent":
        result.append(self.E.adiabaticExponent(self._create_linscalar(prop.text)))
      elif tagname == "refPressure":
        result.append(self.E.refPressure(self._create_linscalar(prop.text)))

    parent.append(result)

  def _convert_testmat(self, testmat, parent):
    result = self.E.testmat()
    
    for prop in testmat:
      if isinstance(prop, etree._Comment):
        parent.append(etree.Comment(prop.text))
        continue
      
      tagname = _clean_ns(prop.tag)
      if tagname == "alpha":
        result.append(self.E.alpha(prop.text))
      elif tagname == "beta":
        result.append(self.E.beta(prop.text))
    
    parent.append(result)

  def _convert_heat(self, heat, parent):
    result = self.E.heatConduction()
    
    for prop in heat:
      if isinstance(prop, etree._Comment):
        result.append(etree.Comment(prop.text))
        continue
      
      tagname = _clean_ns(prop.tag)
      if tagname == "density":
        self._convert_density(prop, result)
      elif tagname == "heatCapacity":
        c = self.E.heatCapacity()
        
        cap_iso = prop.findtext("./linear/isotropic", namespaces=self.nsmap_old)
        if cap_iso is not None:
          c.append(self._create_linscalar(cap_iso))
          
        cap_nl = prop.find("./nonlinear/isotropic", namespaces=self.nsmap_old)
        if cap_nl is not None:
          new_nl = self.E.nonlinear()
          self._convert_nonlindata(new_nl, cap_nl)
          c.append(new_nl)
        
        result.append(c)
      elif tagname == "heatConductivity":
        c = self.E.heatConductivity()
        self._convert_nonlintensor(prop, c)
        result.append(c)
      elif tagname == "refTemperature":
        result.append(self.E.refTemperature(self._create_linscalar(prop.text)))
    
    parent.append(result)

  def _convert_magnetic(self, mag, parent):
    result = self.E.magnetic()
    
    for prop in mag:
      if isinstance(prop, etree._Comment):
        result.append(etree.Comment(prop.text))
        continue
      
      tagname = _clean_ns(prop.tag)
      if tagname == "electricConductivity":
        self._convert_elecconductivity(prop, result)
      elif tagname == "magneticPermeability":
        self._convert_magperm(prop, result)
      elif tagname == "magnetoStrictionTensor_h_mag":
        result.append(self._convert_tensor("magnetoStrictionTensor_h_mag", prop))
      elif (tagname == "DT_ScalarTypeNonNeg" or
            tagname == "hystModel" or
            tagname == "coreLoss"):
        result.append(self._deepcopy(prop))
        
    parent.append(result)

  def _convert_magperm(self, perm, parent):
    result = self.E.permeability()

    lin = perm.find("./linear", namespaces=self.nsmap_old)
    if lin is not None:
      newlin = self.E.linear()
      
      iso = lin.findtext("./isotropic", namespaces=self.nsmap_old)
      if iso is not None:
        newlin.append(self._create_scalar("isotropic", iso))
      
      ortho = lin.find("./orthotropic", namespaces=self.nsmap_old)
      if ortho is not None:
        v1 = ortho.findtext("./permeability_1", namespaces=self.nsmap_old)
        v2 = ortho.findtext("./permeability_2", namespaces=self.nsmap_old)
        v3 = ortho.findtext("./permeability_3", namespaces=self.nsmap_old)
        newlin.append(self._create_orthotensor(v1, v2, v3))
      
      tensor = lin.find("./tensor", namespaces=self.nsmap_old)
      if tensor is not None:
        newlin.append(self._convert_tensor("tensor", tensor))
      
      result.append(newlin)
    
    nonlin = perm.find("./nonlinear", namespaces=self.nsmap_old)
    if nonlin is not None:
      new_nl = self.E.nonlinear()
      
      iso = nonlin.find("./isotropic", namespaces=self.nsmap_old)
      if iso is not None:
        new_iso = self.E.isotropic()
        
        self._convert_nonlindata(new_iso, iso)
        nu = iso.findtext("nuExpr", namespaces=self.nsmap_old)
        if nu is not None:
          new_iso.append(self.E.nuExpr(nu))
        nu_prime = iso.findtext("nuDerivExpr", namespaces=self.nsmap_old)
        if nu_prime is not None:
          new_iso.append(self.E.nuDerivExpr(nu_prime))
        
        new_nl.append(new_iso)
        
      aniso = nonlin.find("anisotropic", namespaces=self.nsmap_old)
      if aniso is not None:
        new_aniso = self.E.anisotropic()
        
        for data in aniso:
          if isinstance(data, etree._Comment):
            result.append(etree.Comment(data.text))
            continue
          
          new_data = self.E.data()
          self._convert_nonlindata(new_data, data)
          nu = data.findtext("nuExpr", namespaces=self.nsmap_old)
          if nu is not None:
            new_data.append(self.E.nuExpr(nu))
          nu_prime = data.findtext("nuDerivExpr", namespaces=self.nsmap_old)
          if nu_prime is not None:
            new_data.append(self.E.nuDerivExpr(nu_prime))
          new_data.append(self.E.angle(data.findtext("angle", namespaces=self.nsmap_old)))
          new_data.append(self.E.zScaling(data.findtext("zScaling", namespaces=self.nsmap_old)))
          
          new_aniso.append(new_data)
          
        new_nl.append(new_aniso)
        
      result.append(new_nl)
        
    parent.append(result)

  def _convert_magstrict(self, magstrict, parent):
    result = self.E.magnetoStrictive()
    
    for prop in magstrict:
      if isinstance(prop, etree._Comment):
        result.append(etree.Comment(prop.text))
        continue
      
      if _clean_ns(prop.tag) == "magnetoStrictionTensor_h":
        result.append(self._convert_tensor("magnetoStrictionTensor_h", prop))
    
    parent.append(result)

  def _convert_piezo(self, piezo, parent):
    result = self.E.piezo()
    
    for prop in piezo:
      if isinstance(prop, etree._Comment):
        result.append(etree.Comment(prop.text))
        continue
      
      tagname = _clean_ns(prop.tag)
      if tagname == "piezoCouplingTensor":
        cpl = self._ensure_elem(result, "piezoCoupling")
        cpl.append(self.E.linear(self._convert_tensor("tensor", prop)))
      elif tagname == "piezoCouplingCoefficient":
        cpl = self._ensure_elem(result, "piezoCoupling")
        self._convert_nonlin(prop, cpl)
      elif tagname == "piezoMicroData":
        result.append(self._deepcopy(prop))
    
    parent.append(result)
  
  def _convert_pyro(self, pyro, parent):
    result = self.E.pyroElectric()
    
    for prop in pyro:
      if isinstance(prop, etree._Comment):
        result.append(etree.Comment(prop.text))
        continue
      
      if _clean_ns(prop.tag) == "pyrocoefficient":
        p = self.E.pyroCoefficient()
        self._convert_tensorprop(prop, p)
        result.append(p)
    
    parent.append(result)
  
  def _convert_thermoelast(self, thermo, parent):
    result = self.E.thermoElastic()
    
    for prop in thermo:
      if isinstance(prop, etree._Comment):
        result.append(etree.Comment(prop.text))
        continue
      
      if _clean_ns(prop.tag) == "thermalExpansion":
        ex = self.E.thermalExpansion()
        self._convert_tensorprop(prop, ex)
        result.append(ex)
        
    parent.append(result)
  
    
if __name__ == "__main__":
  if sys.version_info < (3, 6):
    print("This script requires Python version >= 3.6", file=sys.stderr)
    sys.exit(1)
  
  parser = argparse.ArgumentParser()
  parser.add_argument("input", help="Input XML file using old material schema")
  parser.add_argument("output", help="Output XML file using old material schema")
  
  args = parser.parse_args()
    
  if args.output is None or len(args.output) == 0:
    args.output = args.input
  
  converter = mat_converter()
  ret = converter.convert_mat_xml(args.input, args.output)
  sys.exit(ret)
