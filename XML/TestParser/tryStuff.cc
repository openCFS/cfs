

  // =============================
  //   tryStuff (do experiments)
  // =============================
  void SPHandler::tryStuff() {


    // We try to insert an <output> element, if there is none.
    // The latter may happen, since its occurence is optional.

    bool outputExists = false;

    // ****************************************
    //  STEP1: Check for existance of <output>
    // ****************************************

    // We first determine all <output> elements in the document
    DOMNodeList *outElems = rootelem_->getElementsByTagName( C2X("output") );
    if ( outElems->getLength() == 0 ) {
      outputExists = false;
    }

    else {

      // If there are output elements, we must check, if they are
      // the correct ones:

      // We first locate the <cfsSimulation> element
      DOMNode *cfsSimulation = rootelem_;

      // Now check for each <output> element, if <cfsSimulation> is
      // its parent
      DOMNode *parent;
      for ( unsigned int i = 0; i < outElems->getLength(); i++ ) {
        std::cerr << "Testing <output> node " << i << "\n";
        parent = outElems->item(i)->getParentNode();
        std::cerr << "Parent is <"
                  << X2C( parent->getNodeName() )
                  << ">\n";
        if ( parent->isSameNode(cfsSimulation) ) {
          outputExists = true;
          break;
        }
      }
    }

    std::cerr << "\n --> Status of <output> element: ";
    if ( outputExists ) {
      std::cerr << "exists\n";
    }
    else {
      std::cerr << "does not exist\n";
    }


    // **********************************
    //  STEP2: Generate <output> element
    // **********************************
    if ( outputExists == false ) {
      std::cerr << " --> Generating <output> element\n";

      DOMDocument *doc = parser_->getDocument();
      DOMElement *output = doc->createElementNS(
       C2X( "http://www.cfs++.org" ), C2X("output") );

      // DOMElement *output = doc->createElement( C2X("output") );

      // Check for attributes
      DOMNamedNodeMap *map = output->getAttributes();
      std::cerr << " --> Element has " << map->getLength() << " attributes\n";

      // Append element to tree at correct level
      std::cerr << " --> Appending <output> element to DOM tree\n";
      rootelem_->appendChild( output );

      // Get root grammar
      Grammar* rootGrammar = parser_->getRootGrammar();
      if ( rootGrammar->getGrammarType() == Grammar::SchemaGrammarType ) {
        std::cerr << " --> Got Schema grammar\n";
      }
      else {
        MYABORT;
      }
      SchemaGrammar *sg = dynamic_cast<SchemaGrammar*>(rootGrammar);

      // Get element decl for <output>
      XMLElementDecl *decl = sg->getElemDecl( 29, C2X("output"), C2X(""), Grammar::TOP_LEVEL_SCOPE );
      if ( decl == NULL ) {
        std::cerr << " --> Zefix!\n";
      }

      RefHash3KeysIdPoolEnumerator<SchemaElementDecl> elemEnum = sg->getElemEnumerator();

      // Find output element
      bool gotEmBuggers = false;
      while( elemEnum.hasMoreElements() == true ) {

        const SchemaElementDecl& curElem = elemEnum.nextElement();
        std::string elemName = X2S( curElem.getFullName() );
        // std::cout << "Name: " << elemName << "\n";
        if ( elemName == "output" ) {
          gotEmBuggers = true;
        }
      }

      if ( gotEmBuggers ) {
        std::cerr << "\n\n --> Got 'em Bugger!\n\n";
      }

      process( sg );

   }

    DOMDocument *doc = parser_->getDocument();
    DOMDocumentType *type = doc->getDoctype();
    if ( type == NULL ) {
      std::cerr << "Unable to obtain document type information\n\n";
    }

  }


void SPHandler::process( SchemaGrammar *sg ) {

  RefHash3KeysIdPoolEnumerator<SchemaElementDecl> elemEnum =
    sg->getElemEnumerator();

  if ( !elemEnum.hasMoreElements() ) {
    std::cout << "\nThe validator has no elements to display\n" << std::endl;
    return;
  }

  while( elemEnum.hasMoreElements() ) {

    const SchemaElementDecl& curElem = elemEnum.nextElement();

    // Name
    std::cout << "Name:\t\t\t" << X2C(curElem.getFullName()) << "\n";

    // Model Type
    std::cout << "Model Type:\t\t";
    switch( curElem.getModelType() ) {
    case SchemaElementDecl::Empty:
      std::cout << "Empty";
      break;
    case SchemaElementDecl::Any:
      std::cout << "Any";
      break;
    case SchemaElementDecl::Mixed_Simple:
      std::cout <<
        "Mixed_Simple";  break;
    case SchemaElementDecl::Mixed_Complex:
      std::cout << "Mixed_Complex";
      break;
    case SchemaElementDecl::Children:
      std::cout << "Children";
      break;
    case SchemaElementDecl::Simple:
      std::cout << "Simple";
      break;

    default:
      std::cout << "Unknown";
      break;
    }

    std::cout << "\n";

    // Create Reason
    std::cout << "Create Reason:\t";
    switch( curElem.getCreateReason() ) {
    case XMLElementDecl::NoReason:
      std::cout << "Empty";            break;
    case XMLElementDecl::Declared:
      std::cout << "Declared";         break;
    case XMLElementDecl::AttList:
      std::cout << "AttList";          break;
    case XMLElementDecl::InContentModel:
      std::cout << "InContentModel";   break;
    case XMLElementDecl::AsRootElem:
      std::cout << "AsRootElem";       break;
    case XMLElementDecl::JustFaultIn:
      std::cout << "JustFaultIn";      break;

    default:
      std::cout << "Unknown";  break;
    }

    std::cout << "\n";

    // Content Spec Node
    processContentSpecNode( curElem.getContentSpec() );

    // Misc Flags
    int mflags = curElem.getMiscFlags();
    if( mflags !=0 ) {
      std::cout << "Misc. Flags:\t";
    }

    if ( mflags & SchemaSymbols::XSD_NILLABLE )
      std::cout << "Nillable ";

    if ( mflags & SchemaSymbols::XSD_ABSTRACT )
      std::cout << "Abstract ";

    if ( mflags & SchemaSymbols::XSD_FIXED )
      std::cout << "Fixed ";

    if( mflags !=0 ) {
      std::cout << "\n";
    }

    // Substitution Name
    SchemaElementDecl* subsGroup = curElem.getSubstitutionGroupElem();
    if( subsGroup ) {
      const XMLCh* uriText = parser_->getURIText(subsGroup->getURI());
      std::cout << "Substitution Name:\t" << X2C(uriText)
                << "," << X2C(subsGroup->getBaseName()) << "\n";
    }

    // Content Model
    const XMLCh* fmtCntModel = curElem.getFormattedContentModel();
    if( fmtCntModel != NULL ) {
      std::cout << "Content Model:\t" << X2C(fmtCntModel) << "\n";
    }

    const ComplexTypeInfo* ctype = curElem.getComplexTypeInfo();
    if( ctype != NULL ) {
      std::cout << "ComplexType:\n";
      std::cout << "\tTypeName:\t" << X2C(ctype->getTypeName()) << "\n";

      ContentSpecNode* cSpecNode = ctype->getContentSpec();
      processContentSpecNode(cSpecNode, true );
    }

    // Datatype
    DatatypeValidator* dtValidator = curElem.getDatatypeValidator();
    processDatatypeValidator( dtValidator );

    // Get an enumerator for this guy's attributes if any
    if ( curElem.hasAttDefs() ) {
      processAttributes( curElem.getAttDefList() );
    }

    std::cout << "--------------------------------------------";
    std::cout << std::endl;

  }

  return;
}


//---------------------------------------------------------------------
//  Prints the Attribute's properties
//---------------------------------------------------------------------
void SPHandler::processAttributes( XMLAttDefList& attList, bool margin ) {

  if ( attList.isEmpty() ) {
    return;
  }

  if ( margin ) {
    std::cout << "\t";
  }

  std::cout << "Attributes:\n";
  while( attList.hasMoreElements() ) {

    // Name
    SchemaAttDef& curAttDef = (SchemaAttDef&)attList.nextElement();
    std::cout << "\tName:\t\t\t" << X2C(curAttDef.getFullName()) << "\n";

    // Type
    std::cout << "\tType:\t\t\t";
    std::cout << X2C(XMLAttDef::getAttTypeString(curAttDef.getType()));
    std::cout << "\n";

    // Default Type
    std::cout << "\tDefault Type:\t";
    std::cout << X2C(XMLAttDef::getDefAttTypeString(curAttDef.getDefaultType()));
    std::cout << "\n";

    // Value
    if( curAttDef.getValue() ) {
      std::cout << "\tValue:\t\t\t";
      std::cout << X2C(curAttDef.getValue());
      std::cout << "\n";
    }

    // Enum. values
    if( curAttDef.getEnumeration() ) {
      std::cout << "\tEnumeration:\t";
      std::cout << X2C(curAttDef.getEnumeration());
      std::cout << "\n";
    }

    const DatatypeValidator* dv = curAttDef.getDatatypeValidator();
    processDatatypeValidator( dv, true );

    std::cout << "\n";
  }
}

void SPHandler::processDatatypeValidator( const DatatypeValidator* dtValidator,
                                          bool margin ) {

  if( !dtValidator ) {
    return;
  }

  if( margin ) {
    std::cout << "\t";
  }

  std::cout << "Base Datatype:\t\t";
  switch( dtValidator->getType() ) {
  case DatatypeValidator::String:         std::cout << "string";      break;
  case DatatypeValidator::AnyURI:         std::cout << "AnyURI";      break;
  case DatatypeValidator::QName:          std::cout << "QName";       break;
  case DatatypeValidator::Name:           std::cout << "Name";        break;
  case DatatypeValidator::NCName:         std::cout << "NCName";      break;
  case DatatypeValidator::bool:        std::cout << "bool";     break;
  case DatatypeValidator::Float:          std::cout << "Float";       break;
  case DatatypeValidator::Double:         std::cout << "Double";      break;
  case DatatypeValidator::Decimal:        std::cout << "Decimal";     break;
  case DatatypeValidator::HexBinary:      std::cout << "HexBinary";   break;
  case DatatypeValidator::Base64Binary:   std::cout << "Base64Binary";break;
  case DatatypeValidator::Duration:       std::cout << "Duration";    break;
  case DatatypeValidator::DateTime:       std::cout << "DateTime";    break;
  case DatatypeValidator::Date:           std::cout << "Date";        break;
  case DatatypeValidator::Time:           std::cout << "Time";        break;
  case DatatypeValidator::MonthDay:       std::cout << "MonthDay";    break;
  case DatatypeValidator::YearMonth:      std::cout << "YearMonth";   break;
  case DatatypeValidator::Year:           std::cout << "Year";        break;
  case DatatypeValidator::Month:          std::cout << "Month";       break;
  case DatatypeValidator::Day:            std::cout << "Day";         break;
  case DatatypeValidator::ID:             std::cout << "ID";          break;
  case DatatypeValidator::IDREF:          std::cout << "IDREF";       break;
  case DatatypeValidator::ENTITY:         std::cout << "ENTITY";      break;
  case DatatypeValidator::NOTATION:       std::cout << "NOTATION";    break;
  case DatatypeValidator::List:           std::cout << "List";        break;
  case DatatypeValidator::Union:          std::cout << "Union";       break;
  case DatatypeValidator::AnySimpleType:  std::cout << "AnySimpleType"; break;
  }

  std::cout << "\n";

  // Facets
  RefHashTableOf<KVStringPair>* facets = dtValidator->getFacets();
  if( facets ) {
    RefHashTableOfEnumerator<KVStringPair> enumFacets(facets);
    if( enumFacets.hasMoreElements() ) {
      std::cout << "Facets:\t\t\n";
    }

    while(enumFacets.hasMoreElements()) {
      // Element's properties
      const KVStringPair& curPair = enumFacets.nextElement();
      std::cout << "\t" << X2C( curPair.getKey() )    << "="
                << X2C( curPair.getValue() )  << "\n";
    }
  }

  // Enumerations
  RefVectorOf<XMLCh>* enums = (RefVectorOf<XMLCh>*) dtValidator->getEnumString();
  if (enums) {
    std::cout << "Enumeration:\t\t\n";

    int enumLength = enums->size();
    for ( int i = 0; i < enumLength; i++) {
      std::cout << "\t" << X2C( enums->elementAt(i)) << "\n";
    }
  }
}

void SPHandler::processContentSpecNode( const ContentSpecNode* cSpecNode,
                                        bool margin ) {
  if( !cSpecNode ) {
    return;
  }

  if( margin ) {
    std::cout << "\t";
  }

  std::cout << "ContentType:\t";
  switch( cSpecNode->getType() ) {
  case ContentSpecNode::Leaf:             std::cout << "Leaf";           break;
  case ContentSpecNode::ZeroOrOne:        std::cout << "ZeroOrOne";      break;
  case ContentSpecNode::ZeroOrMore:       std::cout << "ZeroOrMore";     break;
  case ContentSpecNode::OneOrMore:        std::cout << "OneOrMore";      break;
  case ContentSpecNode::Choice:           std::cout << "Choice";         break;
  case ContentSpecNode::Sequence:         std::cout << "Sequence";       break;
  case ContentSpecNode::All:              std::cout << "All";            break;
  case ContentSpecNode::Any:              std::cout << "Any";            break;
  case ContentSpecNode::Any_Other:        std::cout << "Any_Other";      break;
  case ContentSpecNode::Any_NS:           std::cout << "Any_NS";         break;
  case ContentSpecNode::Any_Lax:          std::cout << "Any_Lax";        break;
  case ContentSpecNode::Any_Other_Lax:    std::cout << "Any_Other_Lax";  break;
  case ContentSpecNode::Any_NS_Lax:       std::cout << "Any_NS_Lax";     break;
  case ContentSpecNode::Any_Skip:         std::cout << "Any_Skip";       break;
  case ContentSpecNode::Any_Other_Skip:   std::cout << "Any_Other_Skip"; break;
  case ContentSpecNode::Any_NS_Skip:      std::cout << "Any_NS_Skip";    break;
  case ContentSpecNode::UnknownType:      std::cout << "UnknownType";    break;
  }
  std::cout << "\n";
}
