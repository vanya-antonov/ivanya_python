
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

class BaseAlchemy:
    """Basic parent class to connect to a database using SQLAlchemy.

    Attributes:
     - engine - SQLAlchemy engine
     - s - SQLAlchemy session.
    """

    def __init__(self, db_info):
        """Create SQLAlchemy's engine and session.

        Arguments:
         - db_info - dict with the following keys - host, name, user, pass.
        """
        self.engine = create_engine('sqlite:///:memory:', echo=True)

        Session = sessionmaker(bind=self.engine)
        self.s = Session()

