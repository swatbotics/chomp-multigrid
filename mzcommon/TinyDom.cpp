/*
* Copyright (c) 2008-2014, Matt Zucker
*
* This file is provided under the following "BSD-style" License:
*
* Redistribution and use in source and binary forms, with or
* without modification, are permitted provided that the following
* conditions are met:
*
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
*
* * Redistributions in binary form must reproduce the above
* copyright notice, this list of conditions and the following
* disclaimer in the documentation and/or other materials provided
* with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
* USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/

#include "TinyDom.h"
#include "strutils.h"
#include <expat.h>
#include <vector>
#include <stdexcept>
#include <assert.h>

TinyDomNode::~TinyDomNode() {
  if (_parent) { _parent->removeChild(this); }
}

const TinyDomElement* TinyDomNode::parent() const { return _parent; }
TinyDomElement* TinyDomNode::parent() { return _parent; }

TinyDomNode::TinyDomNode(TinyDomElement* parent): _parent(parent) {
  if (_parent) {
    _parent->addChild(this);
  }
}

void TinyDomNode::reparent(TinyDomElement* parent) {
  if (_parent) { _parent->removeChild(this); }
  _parent = parent;
  if (_parent) { _parent->addChild(this); }
}

//////////////////////////////////////////////////////////////////////

TinyDomElement::TinyDomElement(const std::string& name, TinyDomElement* parent):
  TinyDomNode(parent), _name(name) {}

TinyDomElement::TinyDomElement(const char* name, TinyDomElement* parent):
  TinyDomNode(parent), _name(name) {}

TinyDomElement::TinyDomElement(TinyDomElement* parent):
  TinyDomNode(parent) {}

TinyDomElement::TinyDomElement():
  TinyDomNode(0) {}

TinyDomElement::~TinyDomElement() {
  while (!_children.empty()) {
    delete _children.front();
  }
}

TinyDomNode::NodeType TinyDomElement::type() const { return TypeElement; }

std::string TinyDomElement::attribute(const std::string& key, 
				      const std::string& defaultValue) const {
  for (AttributeList::const_iterator i=_attributes.begin();
       i!=_attributes.end(); ++i) {
    if (i->key == key) {
      return i->data;
    }
  }
  return defaultValue;
}

bool TinyDomElement::haveAttribute(const std::string& key) const {
  for (AttributeList::const_iterator i=_attributes.begin();
       i!=_attributes.end(); ++i) {
    if (i->key == key) {
      return true;
    }
  }
  return false;
}

void TinyDomElement::clearAttribute(const std::string& key) {
  for (AttributeList::iterator i=_attributes.begin();
       i!=_attributes.end(); ++i) {
    if (i->key == key) {
      _attributes.erase(i);
      return;
    }
  }
}
  
void TinyDomElement::setAttribute(const std::string& key, 
				  const std::string& data) {
  TinyDomAttribute newAttr;
  newAttr.key = key;
  newAttr.data = data;
  for (AttributeList::iterator i=_attributes.begin();
       i!=_attributes.end(); ++i) {
    if (i->key == key) {
      *i = newAttr;
      return;
    }
  }
  _attributes.push_back(newAttr);
}

const std::list<TinyDomNode*>& TinyDomElement::children() const {
  return _children;
}

void TinyDomElement::addChild(TinyDomNode* child, 
			      const TinyDomNode* before) {

  if (!child) { return; }
  removeChild(child);

  bool inserted = false;

  if (before) {
    for (ChildList::iterator i=_children.begin(); i!=_children.end(); ++i) {
      if (*i == before) {    
	_children.insert(i, child);
	inserted = true;
	break;
      }
    }
  }

  if (!inserted) { _children.push_back(child); }

  child->_parent = this;

}

void TinyDomElement::removeChild(TinyDomNode* child) {
  if (!child) { return; }
  for (ChildList::iterator i=_children.begin(); i!=_children.end(); ++i) {
    if (*i == child) {
      _children.erase(i);
      child->_parent = 0;
      return;
    }
  }
}

const std::string& TinyDomElement::name() const { return _name; }

void TinyDomElement::setName(const std::string& name) { _name = name; }

const TinyDomElement::AttributeList& TinyDomElement::attributes() const { return _attributes; }

//////////////////////////////////////////////////////////////////////

TinyDomCharacterData::TinyDomCharacterData(TinyDomElement* parent): 
  TinyDomNode(parent) {}

TinyDomCharacterData::TinyDomCharacterData(const std::string& v, 
					   TinyDomElement* parent):
  TinyDomNode(parent), value(v) {}

TinyDomCharacterData::TinyDomCharacterData(const char* v, 
					   TinyDomElement* parent):
  TinyDomNode(parent), value(v) {}

TinyDomCharacterData::~TinyDomCharacterData() {}

TinyDomNode::NodeType TinyDomCharacterData::type() const {
  return TypeCharData;
}

std::string TinyDomCharacterData::trimmedValue() const {
  return trimws(value);
}

  
//////////////////////////////////////////////////////////////////////

struct TinyDomParser {
  TinyDomElement* rootElement;
  std::vector<TinyDomElement*> stack;
  TinyDomParser(): rootElement(0) {}
};

static void commentHandler(void* userData, const XML_Char* data) {
  // no-op
}

static void processingInstructionHandler(void* userData, 
                                         const XML_Char* target, 
                                         const XML_Char* data) {
  // nop
}

static void xmlDeclHandler(void* userData, 
                           const XML_Char* version,
                           const XML_Char* encoding,
                           int standalone) {
  // nop
}

static void startDoctypeDeclHandler(void*, 
                                    const XML_Char*,
                                    const XML_Char*,
                                    const XML_Char*,
                                    int) {
  // nop
}

static void endDoctypeDeclHandler(void*) {
  // nop
}

static void startElementHandler(void *userData,
				const XML_Char *name,
				const XML_Char **atts) {

  TinyDomParser& tdparser = *((TinyDomParser*)userData);

  TinyDomElement* parent = 0;

  if (tdparser.stack.empty()) {
    assert(!tdparser.rootElement);
  } else {
    parent = tdparser.stack.back();
  }

  TinyDomElement* e = new TinyDomElement(name, parent);

  //std::cerr << "creating element with name " << name << " who is child of " << (parent ? parent->name() : "NULL") << "\n";

  for (unsigned int i=0; atts[i]; i+=2) {
    e->setAttribute(atts[i], atts[i+1]);
  }

  if (!tdparser.rootElement) { tdparser.rootElement = e; }

  tdparser.stack.push_back(e);
  
}


static void endElementHandler(void *userData,
			      const XML_Char *name) {

  TinyDomParser& tdparser = *((TinyDomParser*)userData);

  //std::cerr << "closing off element named " << name << "\n";

  assert( !tdparser.stack.empty() );
  
  tdparser.stack.pop_back();

}

static void characterDataHandler(void *userData,
				 const XML_Char *s,
				 int len) {

  TinyDomParser& tdparser = *((TinyDomParser*)userData);

  assert( !tdparser.stack.empty() );

  TinyDomElement* curElement = tdparser.stack.back();

  TinyDomNode* lastChild = 0;

  if (!curElement->children().empty()) {
    lastChild = curElement->children().back();
  }

  if (lastChild && lastChild->type() == TinyDomNode::TypeCharData) {
    // append
    TinyDomCharacterData* chardata = (TinyDomCharacterData*)lastChild;
    chardata->value.append(s, len);
    //std::cerr << "now char data has value '" << chardata->value <<
    //"' and is child of " << curElement->name() << "\n";
  } else {
    new TinyDomCharacterData(std::string(s, len), curElement);
    //std::cerr << "creating char data with value '" <<
    //std::string(s, len) << "' who is child of " << curElement->name() <<
    //"\n";
 }

}


static void defaultHandler(void *userData,
			   const XML_Char *s,
			   int len) {
  for (int i=0; i<len; ++i) {
    if (!isspace(s[i])) {
      throw std::runtime_error("unexpected input: '" + std::string(s, len) + "'");
    }
  }
}







TinyDomElement* TinyDom::parse(std::istream& istr) {

  XML_Parser parser = XML_ParserCreate(0);
  TinyDomParser tdparser;

  XML_SetUserData(parser, &tdparser);
  XML_SetCharacterDataHandler(parser, characterDataHandler);
  XML_SetElementHandler(parser, startElementHandler, endElementHandler);
  XML_SetCommentHandler(parser, commentHandler);
  XML_SetProcessingInstructionHandler(parser, processingInstructionHandler);
  XML_SetXmlDeclHandler(parser, xmlDeclHandler);
  XML_SetDefaultHandler(parser, defaultHandler);
  XML_SetDoctypeDeclHandler(parser, startDoctypeDeclHandler, endDoctypeDeclHandler);

  try {

    while (1) {
      void* buf = XML_GetBuffer(parser, 1024);
      if (!buf) {
	throw std::runtime_error("out of memory!");
      }
      istr.read((char*)buf, 1024);
	  std::streamsize len = istr.gcount();
      if (istr.fail() && !istr.eof()) {
	throw std::runtime_error("failed IO");
      }
      bool isFinal = (istr.eof() || len < 1024);
      if (! XML_ParseBuffer(parser, len, isFinal)) {
	std::ostringstream ostr;
	ostr << "parse error at line " << XML_GetErrorLineNumber(parser)
	     << ", column " << XML_GetErrorColumnNumber(parser) << ": "
	     << XML_ErrorString(XML_GetErrorCode(parser));
	throw std::runtime_error(ostr.str());
				 
      }
      if (isFinal) {
	break;
      }
    }
	
    XML_ParserFree(parser);
	
  } catch (...) {

    //std::cerr << "Got exception: " << e.what() << "\n";

    if (parser) { XML_ParserFree(parser); }
    delete tdparser.rootElement;
    throw;

  }

  return tdparser.rootElement;

}

//////////////////////////////////////////////////////////////////////

std::string TinyDom::escape(const std::string& str) {
  std::string rval = str;
  for (unsigned int pos=0; pos<rval.length(); ++pos) {
    switch (rval[pos]) {
    case '"':
      rval.replace(pos, 1, "&quot;");
      pos += 5;
      break;
    case '<':
      rval.replace(pos, 1, "&lt;");
      pos += 3;
      break;
    case '>':
      rval.replace(pos, 1, "&gt;");
      pos += 3;
      break;
    case '&':
      rval.replace(pos, 1, "&amp;");
      pos += 4;
      break;
    }
  }
  return rval;
}

std::string TinyDom::unescape(const std::string& str) {
  std::string rval = str;
  for (unsigned int pos=0; pos<rval.length(); ++pos) {
    if (rval.substr(pos, 6) == "&quot;") {
      rval.replace(pos, 6, "\"");
    } else if (rval.substr(pos, 4) == "&lt;") {
      rval.replace(pos, 4, "<");
    } else if (rval.substr(pos, 4) == "&gt;") {
      rval.replace(pos, 4, ">");
    } else if (rval.substr(pos, 5) == "&amp;") {
      rval.replace(pos, 5, "&");
    }
  }
  return rval;
}

//////////////////////////////////////////////////////////////////////

static void unparse(const TinyDomNode* node,
		    std::ostream& ostr,
		    unsigned int indent);

static void unparse(const TinyDomElement* element,
		    std::ostream& ostr,
		    unsigned int indent) {

  ostr.flush();

  // print indent
  for (unsigned int i=0; i<indent; ++i) {
    ostr << " ";
  }

  // start of tag
  ostr << "<" << TinyDom::escape(element->name());

  // attributes
  const TinyDomElement::AttributeList& attrs = element->attributes();

  for (TinyDomElement::AttributeList::const_iterator i=attrs.begin(); i!=attrs.end(); ++i) {
    ostr << " " << TinyDom::escape(i->key) << "=\"" << TinyDom::escape(i->data) << "\"";
  }

  TinyDomElement::ChildList::const_iterator j;
  const TinyDomElement::ChildList& children = element->children();
  if (children.empty()) {
    // end of tag
    ostr << "/>\n";
  } else if (children.size() == 1 &&
	     children.front()->type() == TinyDomNode::TypeCharData &&
	     ((TinyDomCharacterData*)children.front())->value.length() < 50) {
    // end of tag
    ostr << ">";
    // first and only child
    ostr << TinyDom::escape(((TinyDomCharacterData*)children.front())->trimmedValue());
    // end tag
    ostr << "</" << TinyDom::escape(element->name()) << ">\n";
  } else {
    // end of start tag
    ostr << ">\n";
    for (j=children.begin(); j!=children.end(); ++j) {
      unparse(*j, ostr, indent+2);
    }
    for (unsigned int i=0; i<indent; ++i) {
      ostr << " ";
    }
    // end of start tag
    ostr << "</" << TinyDom::escape(element->name()) << ">\n";
  }


  ostr.flush();
  
}



void unparse(const TinyDomCharacterData* cdata,
	     std::ostream& ostr,
	     unsigned int indent) {

  ostr.flush();
  
  std::string str = cdata->trimmedValue();
  if (str.empty()) { return; }

  for (unsigned int i=0; i<indent; ++i) {
    ostr << " ";
  }
  ostr << TinyDom::escape(str) << "\n";

  ostr.flush();

}

static void unparse(const TinyDomNode* node,
		    std::ostream& ostr,
		    unsigned int indent) {

  ostr.flush();

  if (node->type() == TinyDomNode::TypeElement) {
    unparse((TinyDomElement*)node, ostr, indent);
  } else {
    unparse((TinyDomCharacterData*)node, ostr, indent);
  }

  ostr.flush();
  
}


void TinyDom::unparse(const TinyDomElement* element,
		      std::ostream& ostr) {

  ::unparse(element, ostr, 0);

}
