import copy
import json
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from datetime import date
from enum import StrEnum
from pathlib import Path
from typing import Any, Self

from Bio import UniProt
from Bio.Seq import Seq

logger = logging.getLogger(Path(__file__).name)


# Found by browsing Tp UniProt taxonomy, starting from root
#   https://www.uniprot.org/taxonomy/160
TP_UNIPROT_ORGANISM_IDS = [
    160,  # Treponema pallidum

    # Subspecies pallidum
    161,  # Treponema pallidum subsp. pallidum (syphilis treponeme)
    243276,  # Treponema pallidum (strain DSM 117211 / Nichols)
    455434,  # Treponema pallidum subsp. pallidum (strain SS14)
    491081,  # Treponema pallidum subsp. pallidum DAL-1
    666714,  # Treponema pallidum subsp. pallidum (strain Chicago)
    686990,  # Treponema pallidum subsp. pallidum str. Mexico A
    1095955,  # Treponema pallidum subsp. pallidum str. Sea 81-4

    # Subspecies pertenue
    168,  # Treponema pallidum subsp. pertenue (Yaws treponeme)
    491078,  # Treponema pallidum subsp. pertenue (strain Samoa D)
    491079,  # Treponema pallidum subsp. pertenue (strain CDC2)
    491080,  # Treponema pallidum subsp. pertenue (strain Gauthier)

    # Subspecies endemicum
    53436,  # Treponema pallidum subsp. endemicum
    1155776,  # Treponema pallidum subsp. endemicum str. Bosnia A
]


ROOT_UNIPROT_URL = "https://www.uniprot.org"


class Database(StrEnum):
    ALPHAFOLDDB = "AlphaFoldDB"
    ARBA = "ARBA"
    BRENDA = "BRENDA"
    CDD = "CDD"
    CH_EBI = "ChEBI"
    DNASU = "DNASU"
    DOI = "DOI"
    DRUG_BANK = "DrugBank"
    EGGNOG = "eggNOG"
    EMBL = "EMBL"
    ENSEMBL_BACTERIA = "EnsemblBacteria"
    ESTHER = "ESTHER"
    EVOLUTIONARY_TRACE = "EvolutionaryTrace"
    FUN_FAM = "FunFam"
    GENE_ID = "GeneID"
    GENE_3D = "Gene3D"
    GO = "GO"
    GOOGLE = "Google"
    HAMAP = "HAMAP"
    HAMAP_RULE = "HAMAP-Rule"
    HOGENOM = "HOGENOM"
    INT_ACT = "IntAct"
    INTERPRO = "InterPro"
    KEGG = "KEGG"
    MDPOSIT = "MDposit"
    MDREPO = "MDRepo"
    MEROPS = "MEROPS"
    NCBI_FAM = "NCBIfam"
    OMA = "OMA"
    ORTHODB = "OrthoDB"
    PANTHER = "PANTHER"
    PATRIC = "PATRIC"
    PDB = "PDB"
    PDB_SUM = "PDBsum"
    PFAM = "Pfam"
    PIR = "PIR"
    PIRNR = "PIRNR"
    PIRSF = "PIRSF"
    PIRSR = "PIRSR"
    PRINTS = "PRINTS"
    PRO = "PRO"
    PROSITE = "PROSITE"
    PROSITE_PRORULE = "PROSITE-ProRule"
    PROTEOMES = "Proteomes"
    PUBMED = "PubMed"
    REFERENCE = "Reference"
    REFSEQ = "RefSeq"
    RHEA = "Rhea"
    RULE_BASE = "RuleBase"
    SABIO_RK = "SABIO-RK"
    SAM = "SAM"
    SFLD = "SFLD"
    SMART = "SMART"
    SMR = "SMR"
    STRING = "STRING"
    SUPFAM = "SUPFAM"
    TCDB = "TCDB"
    UNIPATHWAY = "UniPathway"
    UNIPROTKB = "UniProtKB"


class EntryType(StrEnum):
    REVIEWED = "UniProtKB reviewed (Swiss-Prot)"
    UNREVIEWED = "UniProtKB unreviewed (TrEMBL)"


class ProteinExistence(StrEnum):
    PROTEIN_LEVEL = "1: Evidence at protein level"
    HOMOLOGY = "3: Inferred from homology"
    PREDICTED = "4: Predicted"
    UNCERTAIN = "5: Uncertain"


class NameType(StrEnum):
    RECOMMENDED = "recommendedName"
    SUBMISSION = "submissionNames"
    ALTERNATIVE = "alternativeNames"


class ProteinDescriptionFlag(StrEnum):
    FRAGMENT = "Fragment"
    PRECURSOR = "Precursor"


class CommentType(StrEnum):
    ACTIVITY_REGULATION = "ACTIVITY REGULATION"
    BIOPHYSICOCHEMICAL_PROPERTIES = "BIOPHYSICOCHEMICAL PROPERTIES"
    BIOTECHNOLOGY = "BIOTECHNOLOGY"
    CATALYTIC_ACTIVITY = "CATALYTIC ACTIVITY"
    CAUTION = "CAUTION"
    COFACTOR = "COFACTOR"
    DOMAIN = "DOMAIN"
    FUNCTION = "FUNCTION"
    INDUCTION = "INDUCTION"
    INTERACTION = "INTERACTION"
    MISCELLANEOUS = "MISCELLANEOUS"
    PATHWAY = "PATHWAY"
    PTM = "PTM"
    SEQUENCE_CAUTION = "SEQUENCE CAUTION"
    SIMILARITY = "SIMILARITY"
    SUBCELLULAR_LOCATION = "SUBCELLULAR LOCATION"
    SUBUNIT = "SUBUNIT"


class FeatureType(StrEnum):
    ACTIVE_SITE = "Active site"
    BETA_STRAND = "Beta strand"
    BINDING_SITE = "Binding site"
    CHAIN = "Chain"
    COILED_COIL = "Coiled coil"
    COMPOSITIONAL_BIAS = "Compositional bias"
    DISULFIDE_BOND = "Disulfide bond"
    DNA_BINDING = "DNA binding"
    DOMAIN = "Domain"
    HELIX = "Helix"
    INITIATOR_METHIONINE = "Initiator methionine"
    LIPIDATION = "Lipidation"
    MODIFIED_RESIDUE = "Modified residue"
    MOTIF = "Motif"
    MUTAGENESIS = "Mutagenesis"
    NON_TERMINAL_RESIDUE = "Non-terminal residue"
    REGION = "Region"
    REPEAT = "Repeat"
    SEQUENCE_CONFLICT = "Sequence conflict"
    SEQUENCE_UNCERTAINTY = "Sequence uncertainty"
    SIGNAL = "Signal"
    SITE = "Site"
    TOPOLOGICAL_DOMAIN = "Topological domain"
    TRANSMEMBRANE = "Transmembrane"
    TURN = "Turn"
    ZINC_FINGER = "Zinc finger"


class LocationModifier(StrEnum):
    EXACT = "EXACT"


class SequenceCautionType(StrEnum):
    ERRONEOUS_INITIATION = "Erroneous initiation"
    FRAMESHIFT = "Frameshift"


class DirectionType(StrEnum):
    LEFT_TO_RIGHT = "left-to-right"
    RIGHT_TO_LEFT = "right-to-left"


class CitationType(StrEnum):
    JOURNAL_ARTICLE = "journal article"
    SUBMISSION = "submission"
    THESIS = "thesis"
    UNPUBLISHED_OBSERVATIONS = "unpublished observations"


@dataclass(frozen=True)
class _WebserverSubclass(ABC):
    @classmethod
    @abstractmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs) -> Self:
        pass

    def get_pretty_str(self, **kwargs) -> str:
        """Create a prettified version of this class's attributes.

        Parameters
        ----------
        kwargs
            Optional arguments to specify which attributes to exclude.
            To exclude an attribute, provide the attribute's name as the
            kwarg and set its value to False. These attribute names should be
            namespaced by the name of the class they belong to.

        Examples
        --------
        get_pretty_str(ClassName__custom_attr=False)
            Disables the attribute `custom_attr` from any instances of the
            `ClassName` class.

        Returns
        -------
        str
            Prettified string of this class's attributes.
        """

        pretty: str = ""
        for attribute, value in vars(self).items():
            namespaced_attribute = f"{type(self).__name__}__{attribute}"
            # Handle arguments
            if namespaced_attribute in kwargs:
                arg = kwargs[namespaced_attribute]
                if not isinstance(arg, bool):
                    raise ValueError(
                        f"Unexpected value for kwarg '{namespaced_attribute}': "
                        f"{arg}. See docstring for proper usage."
                    )

                # Skip attributes if user specified as False.
                if not arg:
                    continue

            # Skip empty values
            if value is None:
                continue

            if isinstance(value, _WebserverSubclass):
                value = "\n\t" + value.get_pretty_str(
                    **kwargs
                ).rstrip("\n").replace("\n", "\n\t")
            elif isinstance(value, (list, tuple, set)):
                # Skip empty values
                if len(value) == 0:
                    continue

                new_value = "\n"
                for element in value:
                    if isinstance(element, _WebserverSubclass):
                        new_value += "\t" + element.get_pretty_str(
                            **kwargs
                        ).rstrip("\n").replace("\n", "\n\t")
                    else:
                        new_value += f"\t{element}\n"

                value = new_value.rstrip("\n")
            pretty += self._get_attribute_pretty_str(attribute, value)
        return pretty

    # noinspection method-may-be-static
    def _get_attribute_pretty_str(self, attribute: str, value: Any):
        """Optional method for overriding the default pretty string formatting
        for an attribute.

        Parameters
        ----------
        attribute
            Name of attribute.
        value
            Value of attribute.

        Returns
        -------
        str
            Formatted pretty string for the attribute.
        """

        return f"{attribute.replace("_", " ").title()}: {value}\n"


@dataclass(frozen=True)
class Evidence(_WebserverSubclass):
    evidence_code: str
    source: Database | None
    id: str | None

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            evidence_code=dictionary.pop("evidenceCode"),
            source=Database(dictionary.pop("source")) \
                if "source" in dictionary else None,
            id=dictionary.pop("id", None),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class EntryAudit(_WebserverSubclass):
    first_public_date: date
    last_annotation_update_date: date
    last_sequence_update_date: date
    entry_version: int
    sequence_version: int

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            first_public_date=date.fromisoformat(dictionary.pop("firstPublicDate")),
            last_annotation_update_date=date.fromisoformat(
                dictionary.pop("lastAnnotationUpdateDate")
            ),
            last_sequence_update_date=date.fromisoformat(
                dictionary.pop("lastSequenceUpdateDate")
            ),
            entry_version=int(dictionary.pop("entryVersion")),
            sequence_version=int(dictionary.pop("sequenceVersion")),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class Organism(_WebserverSubclass):
    scientific_name: str
    taxon_id: int
    lineage: tuple[str, ...]
    evidences: tuple[Evidence, ...]
    common_name: str | None

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            scientific_name=dictionary.pop("scientificName"),
            taxon_id=int(dictionary.pop("taxonId")),
            lineage=tuple(dictionary.pop("lineage")),
            evidences=tuple([
                Evidence.from_dict(d)
                for d in dictionary.pop("evidences", [])
            ]),
            common_name=dictionary.pop("commonName", None),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("Organism__evidences", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class Value(_WebserverSubclass):
    value: str
    evidences: tuple[Evidence, ...]
    value_id: str | None

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            value=dictionary.pop("value"),
            evidences=tuple([
                Evidence.from_dict(e)
                for e in dictionary.pop("evidences", [])
            ]),
            value_id=dictionary.pop("id", None)
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("Value__evidences", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class Name(_WebserverSubclass):
    name_type: NameType
    full_name: Value
    short_names: tuple[Value, ...]
    ec_numbers: tuple[Value, ...]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        if "name_type" not in kwargs:
            raise AttributeError("Missing requirement kwarg 'name_type'.")

        inst = cls(
            name_type=kwargs["name_type"],
            full_name=Value.from_dict(dictionary.pop("fullName")),
            short_names=tuple([
                Value.from_dict(d)
                for d in dictionary.pop("shortNames", [])
            ]),
            ec_numbers=tuple([
                Value.from_dict(d)
                for d in dictionary.pop("ecNumbers", [])
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("Name__ec_numbers", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class ProteinDescription(_WebserverSubclass):
    recommended_name: Name | None
    submission_names: tuple[Name, ...]
    alternative_names: tuple[Name, ...]
    flag: ProteinDescriptionFlag | None
    includes: tuple[tuple[Name, ...]]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        includes = dictionary.pop("includes", tuple())
        if len(includes) > 0:
            include_names = []
            for include in includes:
                current_include_names = []
                for name_type, name_datas in include.items():
                    name_type = NameType(name_type)
                    match name_type:
                        case NameType.RECOMMENDED:
                            include_names.append(Name.from_dict(
                                name_datas,
                                name_type=name_type,
                            ))
                        case NameType.SUBMISSION \
                             | NameType.ALTERNATIVE:
                            for name_data in name_datas:
                                include_names.append(
                                    Name.from_dict(
                                        name_data,
                                        name_type=name_type,
                                ))

                if len(current_include_names) > 0:
                    include_names.append(tuple(current_include_names))

            includes = tuple(include_names)

        inst = cls(
            recommended_name=Name.from_dict(
                dictionary.pop("recommendedName"),
                name_type=NameType.RECOMMENDED
            ) if "recommendedName" in dictionary else None,
            alternative_names=tuple([
                Name.from_dict(d, name_type=NameType.ALTERNATIVE)
                for d in dictionary.pop("alternativeNames", {})
            ]),
            submission_names=tuple([
                Name.from_dict(d, name_type=NameType.SUBMISSION)
                for d in dictionary.pop("submissionNames", {})
            ]),
            flag=ProteinDescriptionFlag(dictionary.pop("flag")) \
                if "flag" in dictionary else None,
            includes=includes,
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class Gene(_WebserverSubclass):
    gene_name: Value | None
    ordered_locus_names: tuple[Value, ...]
    orf_names: tuple[Value, ...]
    synonyms: tuple[Value, ...]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            gene_name=Value.from_dict(dictionary.pop("geneName")) \
                if "geneName" in dictionary else None,
            ordered_locus_names=tuple([
                Value.from_dict(d)
                for d in dictionary.pop("orderedLocusNames", [])
            ]),
            orf_names=tuple([
                Value.from_dict(d)
                for d in dictionary.pop("orfNames", [])
            ]),
            synonyms=tuple([
                Value.from_dict(d)
                for d in dictionary.pop("synonyms", [])
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class SubcellularLocation(_WebserverSubclass):
    location: Value
    topology: Value | None
    orientation: Value | None

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            location=Value.from_dict(dictionary.pop("location")),
            topology=Value.from_dict(dictionary.pop("topology")) \
                if "topology" in dictionary else None,
            orientation=Value.from_dict(dictionary.pop("orientation")) \
                if "orientation" in dictionary else None,
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class CrossReference(_WebserverSubclass):
    database: Database
    cross_reference_id: str
    properties: dict[str, str]
    evidences: tuple[Evidence, ...]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            database=Database(dictionary.pop("database")),
            cross_reference_id=dictionary.pop("id"),
            properties={
                item["key"]: item["value"]
                for item in dictionary.pop("properties", [])
            },
            evidences=tuple([
                Evidence.from_dict(d)
                for d in dictionary.pop("evidences", [])
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("CrossReference__evidences", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class Reaction(_WebserverSubclass):
    name: str
    cross_references: tuple[CrossReference, ...]
    evidences: tuple[Evidence, ...]
    ec_number: str | None

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            name=dictionary.pop("name"),
            cross_references=tuple([
                CrossReference.from_dict(d)
                for d in dictionary.pop("reactionCrossReferences", [])
            ]),
            evidences=tuple([
                Evidence.from_dict(d)
                for d in dictionary.pop("evidences", [])
            ]),
            ec_number=dictionary.pop("ecNumber", None),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("Reaction__evidences", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class Cofactor(_WebserverSubclass):
    name: str
    cross_reference: CrossReference
    evidences: tuple[Evidence, ...]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            name=dictionary.pop("name"),
            cross_reference=CrossReference.from_dict(
                dictionary.pop("cofactorCrossReference")),
            evidences=tuple([
                Evidence.from_dict(d)
                for d in dictionary.pop("evidences", [])
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("Cofactor__evidences", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class Velocity(_WebserverSubclass):
    velocity: float
    unit: str
    enzyme: str
    evidences: tuple[Evidence, ...]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            velocity=float(dictionary.pop("velocity")),
            unit=dictionary.pop("unit"),
            enzyme=dictionary.pop("enzyme"),
            evidences=tuple([
                Evidence.from_dict(d)
                for d in dictionary.pop("evidences", [])
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("Velocity__evidences", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class MichaelisConstant(_WebserverSubclass):
    constant: float
    unit: str
    substrate: str
    evidences: tuple[Evidence, ...]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            constant=float(dictionary.pop("constant")),
            unit=dictionary.pop("unit"),
            substrate=dictionary.pop("substrate"),
            evidences=tuple([
                Evidence.from_dict(d)
                for d in dictionary.pop("evidences", [])
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("MichaelisConstant__evidences", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class KineticParameters(_WebserverSubclass):
    maximum_velocities: tuple[Velocity, ...]
    michaelis_constants: tuple[MichaelisConstant, ...]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            maximum_velocities=tuple([
                Velocity.from_dict(d)
                for d in dictionary.pop("maximumVelocities", [])
            ]),
            michaelis_constants=tuple([
                MichaelisConstant.from_dict(d)
                for d in dictionary.pop("michaelisConstants", [])
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class TextsField(_WebserverSubclass):
    texts: tuple[Value, ...]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            texts=tuple([
                Value.from_dict(d)
                for d in dictionary.pop("texts", [])
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class SequenceCaution(_WebserverSubclass):
    sequence_caution_type: SequenceCautionType
    sequence: str
    evidences: tuple[Evidence, ...]
    note: str | None

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            sequence_caution_type=SequenceCautionType(
                dictionary.pop("sequenceCautionType")
            ),
            sequence=dictionary.pop("sequence"),
            evidences=tuple([
                Evidence.from_dict(d)
                for d in dictionary.pop("evidences", [])
            ]),
            note=dictionary.pop("note", None),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("SequenceCaution__evidences", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class Interactant(_WebserverSubclass):
    uniprotkb_accession: str
    int_act_id: str
    gene_name: str | None

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            uniprotkb_accession=dictionary.pop("uniProtKBAccession"),
            int_act_id=dictionary.pop("intActId"),
            gene_name=dictionary.pop("geneName", None),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class Interaction(_WebserverSubclass):
    interactant_one: Interactant
    interactant_two: Interactant
    number_of_experiments: int
    organism_differ: bool

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            interactant_one=Interactant.from_dict(
                dictionary.pop("interactantOne")),
            interactant_two=Interactant.from_dict(
                dictionary.pop("interactantTwo")),
            number_of_experiments=int(dictionary.pop("numberOfExperiments")),
            organism_differ=bool(dictionary.pop("organismDiffer")),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class PhysiologicalReaction(_WebserverSubclass):
    direction_type: DirectionType
    reaction_cross_reference: CrossReference
    evidences: tuple[Evidence, ...]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            direction_type=DirectionType(dictionary.pop("directionType")),
            reaction_cross_reference=CrossReference.from_dict(
                dictionary.pop("reactionCrossReference")),
            evidences=tuple([
                Evidence.from_dict(d)
                for d in dictionary.pop("evidences", [])
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("PhysiologicalReaction__evidences", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class Comment(_WebserverSubclass):
    comment_type: CommentType
    texts: tuple[Value, ...]
    subcellular_locations: tuple[SubcellularLocation, ...]
    reaction: Reaction | None
    cofactors: tuple[Cofactor, ...]
    kinetic_parameters: KineticParameters | None
    ph_dependence: TextsField | None
    note: TextsField | None
    sequence_caution: SequenceCaution | None
    interactions: tuple[Interaction, ...]
    physiological_reactions: tuple[PhysiologicalReaction, ...]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)

        comment_type = CommentType(dictionary.pop("commentType"))
        sequence_caution = None
        if comment_type == CommentType.SEQUENCE_CAUTION:
            sequence_caution = SequenceCaution.from_dict({
                key: dictionary.pop(key) for key in [
                    "sequenceCautionType",
                    "sequence",
                    "note",
                    "evidences",
                ] if key in dictionary
            })

        inst = cls(
            comment_type=comment_type,
            texts=tuple([
                Value.from_dict(d) for d in dictionary.pop("texts", [])
            ]),
            subcellular_locations=tuple([
                SubcellularLocation.from_dict(d)
                for d in dictionary.pop("subcellularLocations", [])
            ]),
            reaction=Reaction.from_dict(dictionary.pop("reaction")) \
                if "reaction" in dictionary else None,
            cofactors=tuple([
                Cofactor.from_dict(d)
                for d in dictionary.pop("cofactors", [])
            ]),
            kinetic_parameters=KineticParameters.from_dict(dictionary.pop(
                "kineticParameters")) if "kineticParameters" in dictionary else None,
            ph_dependence=TextsField.from_dict(
                dictionary.pop("phDependence")) if "phDependence" in dictionary else None,
            note=TextsField.from_dict(dictionary.pop("note")) \
                if "note" in dictionary else None,
            sequence_caution=sequence_caution,
            interactions=tuple([
                Interaction.from_dict(d)
                for d in dictionary.pop("interactions", [])
            ]),
            physiological_reactions=tuple([
                PhysiologicalReaction.from_dict(d)
                for d in dictionary.pop("physiologicalReactions", [])
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class LocationValue(_WebserverSubclass):
    value: int
    modifier: LocationModifier

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            value=dictionary.pop("value"),
            modifier=LocationModifier(dictionary.pop("modifier")),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class Location(_WebserverSubclass):
    start: LocationValue
    end: LocationValue

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            start=LocationValue.from_dict(dictionary.pop("start")),
            end=LocationValue.from_dict(dictionary.pop("end")),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class Ligand(_WebserverSubclass):
    ligand_name: str
    ligand_id: str | None
    label: str | None
    note: str | None

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            ligand_name=dictionary.pop("name"),
            ligand_id=dictionary.pop("id", None),
            label=dictionary.pop("label", None),
            note=dictionary.pop("note", None),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class AlternativeSequence(_WebserverSubclass):
    original_sequence: Seq
    alternative_sequences: tuple[Seq]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        if len(dictionary) == 0:
            return None

        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            original_sequence=Seq(dictionary.pop("originalSequence")),
            alternative_sequences=tuple([
                Seq(s) for s in dictionary.pop("alternativeSequences")
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class Feature(_WebserverSubclass):
    feature_type: FeatureType
    location: Location
    description: str
    evidences: tuple[Evidence, ...]
    feature_id: str | None
    cross_references: tuple[CrossReference, ...]
    ligand: Ligand | None
    alternative_sequence: AlternativeSequence | None

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            feature_type=FeatureType(dictionary.pop("type")),
            location=Location.from_dict(dictionary.pop("location")),
            description=dictionary.pop("description"),
            evidences=tuple([
                Evidence.from_dict(d)
                for d in dictionary.pop("evidences", [])
            ]),
            feature_id=dictionary.pop("featureId", None),
            cross_references=tuple([
                CrossReference.from_dict(d)
                for d in dictionary.pop("featureCrossReferences", [])
            ]),
            ligand=Ligand.from_dict(dictionary.pop("ligand")) \
                if "ligand" in dictionary else None,
            alternative_sequence=AlternativeSequence.from_dict(dictionary.pop(
                "alternativeSequence")) if "alternativeSequence" in dictionary else None,
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("Feature__evidences", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class Keyword(_WebserverSubclass):
    name: str
    category: str
    keyword_id: str
    evidences: tuple[Evidence, ...]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            name=dictionary.pop("name"),
            category=dictionary.pop("category"),
            keyword_id=dictionary.pop("id"),
            evidences=tuple([
                Evidence.from_dict(d)
                for d in dictionary.pop("evidences", [])
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("Keyword__evidences", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class Citation(_WebserverSubclass):
    citation_id: str
    citation_type: CitationType
    authors: tuple[str]
    citation_cross_references: tuple[CrossReference, ...]
    title: str | None
    publication_date: str
    journal: str | None
    first_page: str | None
    last_page: str | None
    volume: int | None
    submission_database: str | None
    institute: str | None
    address: str | None

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            citation_id=dictionary.pop("id"),
            citation_type=CitationType(dictionary.pop("citationType")),
            authors=tuple(dictionary.pop("authors")),
            citation_cross_references=tuple([
                CrossReference.from_dict(d)
                for d in dictionary.pop("citationCrossReferences", [])
            ]),
            title=dictionary.pop("title", None),
            publication_date=dictionary.pop("publicationDate"),
            journal=dictionary.pop("journal", None),
            first_page=dictionary.pop("firstPage", None),
            last_page=dictionary.pop("lastPage", None),
            volume=int(dictionary.pop("volume")) \
                if "volume" in dictionary else None,
            submission_database=dictionary.pop("submissionDatabase", None),
            institute=dictionary.pop("institute", None),
            address=dictionary.pop("address", None),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class ReferenceComment(_WebserverSubclass):
    reference_type: str
    value: str
    evidences: tuple[Evidence, ...]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            reference_type=dictionary.pop("type"),
            value=dictionary.pop("value"),
            evidences=tuple([
                Evidence.from_dict(d)
                for d in dictionary.pop("evidences", [])
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("ReferenceComment__evidences", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class Reference(_WebserverSubclass):
    reference_number: int
    citation: Citation
    reference_positions: tuple[str, ...]
    comments: tuple[ReferenceComment, ...]
    evidences: tuple[Evidence, ...]

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            reference_number=int(dictionary.pop("referenceNumber")),
            citation=Citation.from_dict(dictionary.pop("citation")),
            reference_positions=tuple(dictionary.pop("referencePositions")),
            comments=tuple([
                ReferenceComment.from_dict(p)
                for p in dictionary.pop("referenceComments", [])
            ]),
            evidences=tuple([
                Evidence.from_dict(d)
                for d in dictionary.pop("evidences", [])
            ]),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("Reference__evidences", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class Sequence(_WebserverSubclass):
    value: Seq
    length: int
    mol_weight: int
    crc64: str
    md5: str

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            value=Seq(dictionary.pop("value")),
            length=int(dictionary.pop("length")),
            mol_weight=int(dictionary.pop("molWeight")),
            crc64=dictionary.pop("crc64"),
            md5=dictionary.pop("md5"),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("Sequence__crc64", False)
        kwargs.setdefault("Sequence__md5", False)
        return super().get_pretty_str(**kwargs)


@dataclass(frozen=True)
class ExtraAttributes(_WebserverSubclass):
    count_by_comment_type: dict[CommentType, int]
    count_by_feature_type: dict[FeatureType, int]
    uniparc_id: str

    @classmethod
    def from_dict(cls, dictionary: dict[str, Any], **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            count_by_comment_type={
                CommentType(k): v
                for k, v in dictionary.pop("countByCommentType",{}).items()
            },
            count_by_feature_type={
                FeatureType(k): v
                for k, v in dictionary.pop("countByFeatureType", {}).items()
            },
            uniparc_id=dictionary.pop("uniParcId"),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst


@dataclass(frozen=True)
class Protein(_WebserverSubclass):
    entry_type: EntryType
    primary_accession: str
    secondary_accessions: tuple[str]
    uniprotkb_id: str
    annotation_score: float
    protein_existence: ProteinExistence
    entry_audit: EntryAudit
    organism: Organism
    protein_description: ProteinDescription
    genes: tuple[Gene, ...]
    comments: tuple[Comment, ...]
    features: tuple[Feature, ...]
    keywords: tuple[Keyword, ...]
    references: tuple[Reference, ...]
    uniprotkb_cross_references: tuple[CrossReference, ...]
    sequence: Sequence
    extra_attributes: ExtraAttributes

    @property
    def alternative_names(self) -> list[str]:
        """All alternative names for this protein.

        Includes `gene.gene_name` and all `gene.synonyms`.
        """

        alternative_names = set()
        for gene in self.genes:
            if gene.gene_name:
                alternative_names.add(gene.gene_name.value)
            for synonym in gene.synonyms:
                alternative_names.add(synonym.value)

        return sorted(alternative_names)

    @property
    def aa_seq(self) -> Seq:
        """Wrapper for amino acid sequence found at `self.sequence.value`."""

        return self.sequence.value

    @property
    def uniprotkb_url(self):
        return f"{ROOT_UNIPROT_URL}/uniprotkb/{self.primary_accession}"

    @classmethod
    def from_dict(cls, dictionary: dict, **kwargs):
        dictionary = copy.deepcopy(dictionary)
        inst = cls(
            entry_type=EntryType(dictionary.pop("entryType")),
            primary_accession=dictionary.pop("primaryAccession"),
            secondary_accessions=tuple(dictionary.pop("secondaryAccessions", [])),
            uniprotkb_id=dictionary.pop("uniProtkbId"),
            entry_audit=EntryAudit.from_dict(dictionary.pop("entryAudit")),
            annotation_score=dictionary.pop("annotationScore"),
            organism=Organism.from_dict(dictionary.pop("organism")),
            protein_existence=ProteinExistence(
                dictionary.pop("proteinExistence")
            ),
            protein_description=ProteinDescription.from_dict(
                dictionary.pop("proteinDescription")),
            genes=tuple([
                Gene.from_dict(d) for d in dictionary.pop("genes", [])
            ]),
            comments=tuple([
                Comment.from_dict(d) for d in dictionary.pop("comments", [])
            ]),
            features=tuple([
                Feature.from_dict(d) for d in dictionary.pop("features", [])
            ]),
            keywords=tuple([
                Keyword.from_dict(d) for d in dictionary.pop("keywords", [])
            ]),
            references=tuple([
                Reference.from_dict(d)
                for d in dictionary.pop("references", [])
            ]),
            uniprotkb_cross_references=tuple([
                CrossReference.from_dict(d)
                for d in dictionary.pop("uniProtKBCrossReferences", [])
            ]),
            sequence=Sequence.from_dict(dictionary.pop("sequence")),
            extra_attributes=ExtraAttributes.from_dict(
                dictionary.pop("extraAttributes")),
        )
        if len(dictionary) > 0:
            raise ValueError(
                f"Remaining data:\n{json.dumps(dictionary, indent=2)}"
            )

        return inst

    def get_pretty_str(self, **kwargs) -> str:
        kwargs.setdefault("Protein__entry_audit", False)
        kwargs.setdefault("Protein__organism", True)
        kwargs.setdefault("Protein__protein_description", True)
        kwargs.setdefault("Protein__genes", True)
        kwargs.setdefault("Protein__comments", False)
        kwargs.setdefault("Protein__features", False)
        kwargs.setdefault("Protein__keywords", False)
        kwargs.setdefault("Protein__references", False)
        kwargs.setdefault("Protein__uniprotkb_cross_references", False)
        kwargs.setdefault("Protein__sequence", True)
        kwargs.setdefault("Protein__extra_attributes", False)
        return super().get_pretty_str(**kwargs)


def get_all_uniprot_tp_proteins() -> list[Protein]:
    uniprot_results = []
    for organism_id in TP_UNIPROT_ORGANISM_IDS:
        logger.debug(f"Searching UniProt for organism {organism_id}...")
        for result in UniProt.search(f"organism_id:{organism_id}"):
            uniprot_results.append(result)

    proteins: list[Protein] = []
    for result in uniprot_results:
        proteins.append(Protein.from_dict(result))

    return proteins
