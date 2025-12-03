from sqlalchemy import Column, Integer, String, Text, create_engine
from sqlalchemy.orm import sessionmaker, declarative_base
from pathlib import Path

# Set up base for SQLAlchemy models
Base = declarative_base()

# Get the directory where this module is located
MODULE_DIR = Path(__file__).parent.parent
DB_DIR = MODULE_DIR / "db"
DB_PATH = DB_DIR / "smart_chemist.db"


class AnnotatedPattern(Base):
    __tablename__ = "annotated_pattern"
    id = Column(Integer, primary_key=True, autoincrement=True)
    trivial_name = Column(String(100), nullable=False, unique=False)
    smarts = Column(Text, nullable=False)
    group = Column(String(50), nullable=True)
    index_file = Column(Integer, nullable=True)
    heavy_atoms = Column(Integer, nullable=True, default=0)
    num_rings = Column(Integer, nullable=True, default=0)
    n_nitrogens = Column(Integer, nullable=True, default=0)
    n_oxygen = Column(Integer, nullable=True, default=0)
    n_sulfur = Column(Integer, nullable=True, default=0)
    n_carbon = Column(Integer, nullable=True, default=0)
    n_halogens = Column(Integer, nullable=True, default=0)
    n_phosphor = Column(Integer, nullable=True, default=0)
    n_other_atom = Column(Integer, nullable=True, default=0)
    hierarchy = Column(String(100), nullable=True)


# Create engine and session
engine = create_engine(f"sqlite:///{DB_PATH}")
Session = sessionmaker(bind=engine)


def get_session():
    """Get a new database session."""
    return Session()


# Export all models and utilities
__all__ = ["AnnotatedPattern", "Base", "engine", "Session", "get_session", "DB_PATH"]

